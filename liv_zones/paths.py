# %%
# Edit these paths to point to your data directories.
# Each path should be an acinus directory containing stack subdirectories.
#
# Expected directory structure:
#   <root>/<Sex>/<Diet>/<Liver>/<Lobule>/<Acinus>/
#
# Example:
#   /path/to/your/data/Male/CNT/Liv1/Lobule1/acinus0/
#
# Diet abbreviations used in this study:
#   CNT = control, STV = steatohepatitis, WD = western diet

male_cnt_paths = [
    '/path/to/your/data/Male/CNT/Liv1/Lobule1/acinus0',
    '/path/to/your/data/Male/CNT/Liv1/Lobule1/acinus1',
    # add more paths as needed...
]

male_stv_paths = [
    '/path/to/your/data/Male/STV/Liv1/Lobule1/acinus0',
    '/path/to/your/data/Male/STV/Liv1/Lobule1/acinus1',
    # add more paths as needed...
]

male_wd_paths = [
    '/path/to/your/data/Male/WD/Liv1/Lobule1/acinus0',
    '/path/to/your/data/Male/WD/Liv1/Lobule1/acinus1',
    # add more paths as needed...
]

female_cnt_paths = [
    '/path/to/your/data/Female/CNT/Liv1/Lobule1/acinus0',
    '/path/to/your/data/Female/CNT/Liv1/Lobule1/acinus1',
    # add more paths as needed...
]

female_stv_paths = [
    '/path/to/your/data/Female/STV/Liv1/Lobule1/acinus0',
    '/path/to/your/data/Female/STV/Liv1/Lobule1/acinus1',
    # add more paths as needed...
]

female_wd_paths = [
    '/path/to/your/data/Female/WD/Liv1/Lobule1/acinus0',
    '/path/to/your/data/Female/WD/Liv1/Lobule1/acinus1',
    # add more paths as needed...
]

liver_paths = male_cnt_paths + male_stv_paths + male_wd_paths + female_cnt_paths + female_stv_paths + female_wd_paths
# %%
