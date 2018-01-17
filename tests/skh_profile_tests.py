"""
The purpose of this script is to profile the length of `phydms` tests.

SKH 20170901
"""
import time
from datetime import timedelta
import glob
import subprocess
import pandas as pd

def main():
    df = {"test":[], "time":[]}
    for test in glob.glob("test_*py"):
        df["test"].append(test)
        cmd = ["python3", test]
        start_time = time.time()
        subprocess.check_call(cmd)
        elapsed_time = time.time() - start_time
        df["time"].append(str(timedelta(seconds=elapsed_time)))
    df = pd.DataFrame(df)
    df.to_csv("skh_profile.csv", index=False)

if __name__ == '__main__':
    main()
