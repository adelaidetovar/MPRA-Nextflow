import os
import json
import argparse

def list_of_strings(arg):
    return arg.split(',')

def make_config(libnames, dirs):
    config = {
        "libnames": {}
    }

    umi_len_map = {
        "_RNA": 15,
        "_cDNA": 15,
        "_DNA": 13,
        "_gDNA": 13
    }

    for libname in libnames:
        umi_len = None
        for pattern, length in umi_len_map.items():
            if pattern in libname:
                umi_len = length
                break

        if umi_len is None:
            raise ValueError(f"Check naming convention for sample '{libname}'")

        lib_entry = {
            "umi_len": [umi_len],
            "fwd": {},
            "rev": {}
        }

        entry = 1

        for dir in dirs:
            fwd_path = os.path.join(dir, f"{libname}.r1.fq.gz")
            rev_path = os.path.join(dir, f"{libname}.r2.fq.gz")

            lib_entry["fwd"][str(entry)] = fwd_path
            lib_entry["rev"][str(entry)] = rev_path

            entry += 1

        config["libnames"][libname] = lib_entry

    return config

def main():
    opts = argparse.ArgumentParser()

    opts.add_argument("--libnames", dest="libnames",type=list_of_strings)
    opts.add_argument("--dirs", dest = "dirs",type=list_of_strings)

    args = opts.parse_args()

    config = make_config(args.libnames, args.dirs)

    with open("input-config.json", "w") as outfile:
        json.dump(config, outfile, indent=4)


if __name__ == "__main__":
    main()
