import logging
import subprocess
from typing import List, Dict

import pandas as pd


def get_filtered_dataframe(file: str, cols_to_compare: Dict[str, type]) -> pd.DataFrame:
    return pd.read_csv(
        file,
        sep="\t",
        dtype=cols_to_compare,
        usecols=cols_to_compare.keys(),
    )


def run_cmd(cmd: str, raise_ex: bool = True) -> bool:
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error: {e}")
        if raise_ex:
            raise
        return False
    return True


def get_biosample_names(biosamples_tsv: str) -> List[str]:
    biosample_data = pd.read_csv(biosamples_tsv, sep="\t")
    return biosample_data.iloc[:, 0].tolist()
