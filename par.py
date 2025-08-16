# ============================================================
# Synthetic Time-Series Generation with DeepEcho (PARModel)
# Purpose:
#   - Load and preprocess patient telemetry data.
#   - Train a PARModel per patient.
#   - Sample synthetic sequences per patient and merge results.
#
# Notes:
#   - Preserves original logic and thresholds.
#   - Fixes an original column bug: use 'Datetime_character' (not 'Datetime').
#   - Keeps entity_columns == ['ID'] as intended.
#   - CUDA is auto-detected; set FORCE_CUDA to override.
# ============================================================

# (Optional) installs:
# !pip install deepecho
# !pip install sdv
import matplotlib.pyplot as plt

from __future__ import annotations
import os
import warnings
import pandas as pd

# Third-party
from deepecho import PARModel

# ----------------------------- Config ---------------------------------
DATA_PATH = "..."      # original path
MIN_ROWS_PER_ID = 30                       # keep IDs with >= 30 rows
FILTERS = {
    "Flow_min": 2.0,                       # df['Flow'] >= 2
    "Motor_power_min": 2.5,                # df['Motor_power'] >= 2.5
}
DIVIDE_STORED_SPEED_BY = 1000.0            # Stored_speed / 1000
SEQUENCE_INDEX = "Datetime_character"      # used by DeepEcho
ENTITY_COLUMNS = ["ID"]             # DeepEcho grouping key
NUM_ENTITIES_PER_PATIENT = 50              # synthetic sequences per patient
SEQ_LENGTH = 360                           # synthetic sequence length
EPOCHS_GLOBAL = 1000                       # (kept for parity; not used in per-patient loop)
EPOCHS_PER_PATIENT = 500                   # epochs per patient model
FORCE_CUDA = None                          # None => auto-detect; True/False to override

# Data types required by deepecho (kept as in original)
DATA_TYPES = {
    "Abbott_ID": "categorical",
    "HCT": "categorical",
    "Flow": "continuous",
    "Motor_power": "continuous",
    "Stored_speed": "categorical",
    "Hour": "categorical",
    "Date": "categorical",
}

# ----------------------- Utilities & Helpers ---------------------------
def _cuda_flag() -> bool:
    """Decide CUDA usage. Auto-detect (torch) unless FORCE_CUDA specified."""
    if FORCE_CUDA is not None:
        return bool(FORCE_CUDA)
    try:
        import torch  # noqa: F401
        import torch.cuda
        return torch.cuda.is_available()
    except Exception:
        return False


def load_and_preprocess(path: str) -> pd.DataFrame:
    """
    Load CSV, select columns, filter, handle NA, transform speed, and derive time columns.
    Mirrors original logic; fixes 'Datetime' -> 'Datetime_character' reference.
    """
    # Load
    df = pd.read_csv(path)

    # Keep columns of interest (original selection)
    df = df[["ID", "Datetime", "Motor_power", "Flow", "HCT", "Stored_speed"]]

    # Keep IDs with at least MIN_ROWS_PER_ID rows
    ids_with_min_n = (
        df.groupby("ID").filter(lambda x: len(x) >= MIN_ROWS_PER_ID)["ID"].unique()
    )
    df = df[df.idxmax.isin(ids_with_min_n)].copy()

    # Drop duplicates
    df.drop_duplicates(inplace=True)


    # Fill NA with numeric column means (original behavior)
    # (Pandas: mean() on DF returns numeric columns' means)
    df.fillna(df.mean(numeric_only=True), inplace=True)

    # Remove weird values (original thresholds)
    df = df[df["Flow"] >= FILTERS["Flow_min"]]
    df = df[df["Motor_power"] >= FILTERS["Motor_power_min"]]

    # Scale Stored_speed
    df["Stored_speed"] = df["Stored_speed"].apply(lambda x: x / DIVIDE_STORED_SPEED_BY)

    # Parse datetime
    df[SEQUENCE_INDEX] = pd.to_datetime(df["Datetime_character"])

    # Derive Date/Hour from Datetime_character (original intended behavior)
    df["Date"] = df[SEQUENCE_INDEX].dt.date
    df["Hour"] = df[SEQUENCE_INDEX].dt.hour

    # Lag-difference of Stored_speed (original: global shift, not grouped)
    df = df.sort_values([ "idxmax", SEQUENCE_INDEX ]).reset_index(drop=True)
    df["Stored_speed_lag"] = df["Stored_speed"].shift(1)
    df["Stored_speed"] = df["Stored_speed"] - df["Stored_speed_lag"]
    df.drop(columns=["Stored_speed_lag"], inplace=True)

    return df


class DeepEchoSynthesizer:
    """
    Wrapper for training and sampling PARModel per patient.
    Preserves the original logic:
      - Fit a new PARModel per patient (epochs=EPOCHS_PER_PATIENT, cuda=auto).
      - Sample `NUM_ENTITIES_PER_PATIENT` sequences of length SEQ_LENGTH per patient.
      - Offset idxmaxs per patient block to remain unique in the synthetic output.
    """

    def __init__(
        self,
        data_types: dict,
        entity_columns: list[str],
        sequence_index: str,
        epochs_per_patient: int = EPOCHS_PER_PATIENT,
        num_entities: int = NUM_ENTITIES_PER_PATIENT,
        seq_length: int = SEQ_LENGTH,
    ):
        self.data_types = data_types
        self.entity_columns = entity_columns
        self.sequence_index = sequence_index
        self.epochs_per_patient = epochs_per_patient
        self.num_entities = num_entities
        self.seq_length = seq_length
        self.use_cuda = _cuda_flag()

    def fit_and_sample_per_patient(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        For each idxmax, train a PARModel and sample synthetic data.
        Returns a single concatenated DataFrame of all synthetic rows.
        """
        unique_ids = df.idxmax.unique()
        data = []
        block_index = 0

        for pid in unique_ids:
            patient = df[df.idxmax == pid].reset_index(drop=True)

            # Train per-patient model
            model = PARModel(epochs=self.epochs_per_patient, cuda=self.use_cuda)
            model.fit(
                data=patient,
                entity_columns=self.entity_columns,
                data_types=self.data_types,
                sequence_index=self.sequence_index,
            )

            # Sample synthetic sequences
            synth = model.sample(num_entities=self.num_entities, sequence_length=self.seq_length)

            # Maintain unique idxmax blocks in synthetic output (original pattern)
            synth["idxmax"] = synth["idxmax"] + block_index * self.num_entities + 1
            block_index += 1

            print(f"Trained + sampled for patient idxmax={pid}")
            data.append(synth)

            # free references
            del synth, model

        # Combine
        if not data:
            return pd.DataFrame(columns=df.columns)
        return pd.concat(data, axis=0, ignore_index=True)


def main():
    warnings.filterwarnings("ignore")

    # Load and preprocess
    df = load_and_preprocess(DATA_PATH)

    # Re-assert required typing for DeepEcho (just to be explicit)
    # Note: DeepEcho handles internal encoding of categorical/continuous types.
    #       We provide the mapping via DATA_TYPES.
    synthesizer = DeepEchoSynthesizer(
        data_types=DATA_TYPES,
        entity_columns=ENTITY_COLUMNS,
        sequence_index=SEQUENCE_INDEX,
        epochs_per_patient=EPOCHS_PER_PATIENT,
        num_entities=NUM_ENTITIES_PER_PATIENT,
        seq_length=SEQ_LENGTH,
    )

    # Fit per patient and sample synthetic data
    synth_all = synthesizer.fit_and_sample_per_patient(df)

    # Example: save to disk (optional)
    out_path = "..."
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    synth_all.to_csv(out_path, index=False)
    print(f"Synthetic data saved to: {out_path}")

    # ------------- Optional plotting (as in your original, commented) ---------
    pid_min = synth_all.idxmax.min()
    patient_1 = synth_all[synth_all.idxmax == pid_min].reset_index(drop=True)
    real_patient_1 = df[df.idxmax == df.idxmax.unique()[0]].reset_index(drop=True)
    
    fig, ax = plt.subplots(4, 1, figsize=(16, 10))
    ax[0].set_title(f'Patient {pid_min} Motor Power')
    ax[0].plot(patient_1.Motor_power[:100], label='Synthetic')
    ax[0].plot(real_patient_1.Motor_power[:100], label='Real', alpha=0.5)
    
    ax[1].set_title(f'Patient {pid_min} Flow')
    ax[1].plot(patient_1.Flow[:100], label='Synthetic')
    ax[1].plot(real_patient_1.Flow[:100], label='Real', alpha=0.5)
    
    ax[2].set_title(f'Patient {pid_min} HCT')
    ax[2].plot(patient_1.HCT[:100], label='Synthetic')
    ax[2].plot(real_patient_1.HCT[:100], label='Real', alpha=0.5)
    
    ax[3].set_title(f'Patient {pid_min} Stored_speed (Δ)')
    ax[3].plot(patient_1.Stored_speed[:100], label='Synthetic')
    ax[3].plot(real_patient_1.Stored_speed[:100], label='Real', alpha=0.5)
    
    for a in ax: a.legend()
    plt.tight_layout()
    plt.savefig('simulation_figures/pws_deepecho_example.png')
    plt.close(fig)


if __name__ == "__main__":
    main()
