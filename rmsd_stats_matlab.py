#!/usr/bin/env python3

"""
Created on 07 May. 2024
"""

__author__ = "Nicolas JEANNE"
__copyright__ = "GNU General Public License"
__email__ = "jeanne.n@chu-toulouse.fr"
__version__ = "1.0.0"

import argparse
import logging
import os
import re
import shutil
import sys

import pandas as pd


def create_log(path, level):
    """Create the log as a text file and as a stream.

    :param path: the path of the log.
    :type path: str
    :param level: the level og the log.
    :type level: str
    :return: the logging:
    :rtype: logging
    """

    log_level_dict = {"DEBUG": logging.DEBUG,
                      "INFO": logging.INFO,
                      "WARNING": logging.WARNING,
                      "ERROR": logging.ERROR,
                      "CRITICAL": logging.CRITICAL}

    if level is None:
        log_level = log_level_dict["INFO"]
    else:
        log_level = log_level_dict[level]

    if os.path.exists(path):
        os.remove(path)

    logging.basicConfig(format="%(asctime)s %(levelname)s:\t%(message)s",
                        datefmt="%Y/%m/%d %H:%M:%S",
                        level=log_level,
                        handlers=[logging.FileHandler(path), logging.StreamHandler()])
    return logging


def get_conditions(path):
    """
    Extract the conditions, the paths and the colors.

    :param path: the path to the CSV file.
    :type path: str
    :return: the conditions.
    :rtype: pd.DataFrame
    """
    df = pd.read_csv(path, sep=",", header=None)
    df.columns = ["condition", "path", "color"]
    return df


def extract_rmsd(conditions, index_file, data_out_path):
    """
    Extract the RMSD values of each sample and return the aggregated RMSD values for each frame.

    :param conditions: the conditions dataframe.
    :type conditions: pandas.DataFrame
    :param index_file: the path to the index file.
    :type index_file: str
    :param data_out_path: the path to output data file.
    :type data_out_path: str
    :return: the regrouped data.
    :rtype: pandas.DataFrame
    """

    index_samples = pd.read_excel(index_file, sheet_name=0)
    index_conditions = pd.read_excel(index_file, sheet_name=1)

    data = {"sample": [], "condition": [], "frame": [],  f"RMSD": []}
    pattern = re.compile("RMSD_(.+)_.+\\.csv")
    conditions_to_remove = []
    for _, row_condition in conditions.iterrows():
        by_condition = [fn for fn in os.listdir(row_condition["path"]) if
                        fn.startswith("RMSD") and not "histogram" in fn and fn.endswith(".csv")]
        if len(by_condition) == 0:
            conditions_to_remove.append(row_condition["condition"])
            logging.warning(f"Condition {row_condition['condition']}: no RMSD files, this condition is skipped.")
            continue
        logging.info(f"Retrieving {len(by_condition)} file{'s' if len(by_condition) > 1 else ''} data for condition: "
                     f"{row_condition['condition']}")
        for item in sorted(by_condition):
            match = pattern.search(item)
            if match:
                sample = match.group(1)
            else:
                logging.error(f"\tNo match between the pattern '{pattern.pattern}' and the file name {item}.")
                sys.exit(1)
            df_current = pd.read_csv(os.path.join(row_condition["path"], item), sep=",")
            prefix = f"\t\t- {sample}:"
            logging.info(f"{prefix:<50}{len(df_current['frames'])} frames.")
            current_frames = df_current["frames"].to_list()
            data["sample"] = data["sample"] + (
                        [index_samples[index_samples["sample"] == sample]["index"].values[0]] * len(current_frames))
            try:
                data["condition"] = data["condition"] + (
                        [index_conditions[index_conditions["condition"] == row_condition["condition"]]["index"].values[
                            0]] * len(current_frames))
            except IndexError as ex:
                logging.error(f"\"{row_condition['condition']}\" does not exists in the second tab of the index file: "
                              f"{index_file}")
                sys.exit(1)
            data["frame"] = data["frame"] + current_frames
            data["RMSD"] = data["RMSD"] + df_current["RMSD"].to_list()

    df = pd.DataFrame.from_dict(data)
    df.to_csv(data_out_path, sep=',', index=False)
    logging.info(f"RMSD data by sample and condition written: {data_out_path}")

    destination_index_file = os.path.join(os.path.dirname(data_out_path), os.path.basename(index_file))
    shutil.copyfile(index_file, destination_index_file)
    logging.info(f"Index file copied to: {destination_index_file}")


if __name__ == "__main__":
    descr = f"""
    {os.path.basename(__file__)} v. {__version__}

    Created by {__author__}.
    Contact: {__email__}
    {__copyright__}

    Distributed on an "AS IS" basis without warranties or conditions of any kind, either express or implied.

    From the RMSD data by sample, create a CSV file of all the RMSD values of all the samples recoded for MatLab usage. 
    The RMS files are the RMSD CSV files from the RMS analysis (https://github.com/njeanne/rms). 
    The Molecular Dynamics simulation time must be identical.

    The inputs are:
    
    - a comma separated file without header which first column is the condition, the second column the path of the 
    directory containing the RMS analysis files and the third column the color in hexadecimal format. i.e:

    insertions,tests/inputs/insertions,#fc030b
    WT,tests/inputs/WT,#0303fc
    
    - an XLS file with 2 tabs, one with the samples and the corresponding index, the other with the conditions and the 
    index corresponding. The first tab contain the samples and the corresponding index, the second tab the conditions 
    and the corresponding index. See data/index.xlsx file for more details.

    The output is CSV file of all the samples RMSD data.
    """
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--out", required=True, type=str, help="the path to the output file.")
    parser.add_argument("-i", "--index", required=True, type=str, help="the path to the XLS index file.")
    parser.add_argument("-l", "--log", required=False, type=str,
                        help="the path for the log file. If this option is skipped, the log file is created in the "
                             "output directory.")
    parser.add_argument("--log-level", required=False, type=str,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="set the log level. If the option is skipped, log level is INFO.")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("input", type=str,
                        help="the path to the CSV (comma separated without header) file which first column is the "
                             "condition, the second column the path of the directory containing the RMS analysis "
                             "files and the third column the color.")
    args = parser.parse_args()

    # create output directory if necessary
    out_dir = os.path.dirname(args.out)
    os.makedirs(out_dir, exist_ok=True)
    # create the logger
    if args.log:
        log_path = args.log
    else:
        log_path = os.path.join(out_dir, f"{os.path.splitext(os.path.basename(__file__))[0]}.log")
    create_log(log_path, args.log_level)

    logging.info(f"version: {__version__}")
    logging.info(f"CMD: {' '.join(sys.argv)}")

    data_conditions = get_conditions(args.input)
    extract_rmsd(data_conditions, args.index, args.out)
