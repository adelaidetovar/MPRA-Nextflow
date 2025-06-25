import os
import re
import json
import argparse
from collections import defaultdict

def list_of_strings(arg):
    return arg.split(',')

# change if UMI lengths do not differ by sample type
def get_umi_length(libname):
    umi_len_map = {
        "_RNA": 15,
        "_cDNA": 15,
        "_DNA": 13,
        "_gDNA": 13
    }
    for pattern, length in umi_len_map.items():
        if pattern in libname:
            return length
    return None

def make_config(libnames, dirs, single_end=False):
    config = {
        "libnames": {}
    }

    fq_files = []
    for dir in dirs:
        fq_files.extend([os.path.join(dir, f) for f in os.listdir(dir) if f.endswith(".fq.gz")])

    for libname in libnames:
        umi_len = get_umi_length(libname)
        if not single_end and umi_len is None:
            raise ValueError(f"Check naming convention for sample '{libname}.'")

        fwd_dict = {}
        rev_dict = {} if not single_end else None

        matches = []
        for fq in fq_files:
            fname = os.path.basename(fq)
            if fname.startswith(libname) and re.search(r'_[A-Z]\.r1\.fq\.gz$', fname):
                suffix = re.search(r'_([A-Z])\.r1\.fq\.gz$', fname).group(1)
                fwd = fq
                rev = fq.replace(".r1.fq.gz", ".r2.fq.gz")
                if not single_end and not os.path.exists(rev):
                    raise FileNotFoundError(f"Missing reverse read: {rev}")
                matches.append((suffix, fwd, rev if not single_end else None))

        matches.sort()

        for i, (suffix, fwd, rev) in enumerate(matches, start=1):
            fwd_dict[str(i)] = fwd
            if not single_end:
                rev_dict[str(i)] = rev

        entry = {
            "fwd": fwd_dict
        }

        if not single_end:
            entry["rev"] = rev_dict
            entry["umi_len"] = umi_len

        config["libnames"][libname] = entry

    return config

def main():
    opts = argparse.ArgumentParser()

    opts.add_argument("--libnames", dest="libnames", type=list_of_strings, required=True)
    opts.add_argument("--dirs", dest="dirs", type=list_of_strings, required=True)
    opts.add_argument("--no_umi", action="store_true", help="Use this flag if you're not using UMIs")

    args = opts.parse_args()

    config = make_config(args.libnames, args.dirs, single_end=args.no_umi)

    with open("input-config.json", "w") as outfile:
        json.dump(config, outfile, indent=4)

if __name__ == "__main__":
    main()
