import os
import re
import json
import argparse

def list_of_strings(arg):
    return arg.split(',')

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

def make_config(libnames, dirs, single_end=False, ignore_umi=False):
    config = {
        "libnames": {}
    }

    fq_files = []
    for dir in dirs:
        fq_files.extend([os.path.join(dir, f) for f in os.listdir(dir)
                         if re.search(r"\.r[12]\.(fq|fastq)\.gz$", f)])

    for libname in libnames:
        umi_len = get_umi_length(libname) if not ignore_umi else None
        if not ignore_umi and umi_len is None:
            raise ValueError(f"Check naming convention for sample '{libname}'")

        fwd_dict = {}
        rev_dict = {} if not single_end else None

        matched = []
        for fq in fq_files:
            fname = os.path.basename(fq)
            ext = ".fq.gz" if ".fq.gz" in fname else ".fastq.gz"
            rep_match = re.match(rf"^{re.escape(libname)}_([A-Z])\.r1\.(fq|fastq)\.gz$", fname)
            if rep_match:
                suffix = rep_match.group(1)
                fwd = fq
                rev = fq.replace(".r1", ".r2") if not single_end else None
                if rev and not os.path.exists(rev):
                    raise FileNotFoundError(f"Missing reverse read: {rev}")
                matched.append((suffix, fwd, rev))
            elif re.match(rf"^{re.escape(libname)}\.r1\.(fq|fastq)\.gz$", fname):
                fwd = fq
                rev = fq.replace(".r1", ".r2") if not single_end else None
                if rev and not os.path.exists(rev):
                    raise FileNotFoundError(f"Missing reverse read: {rev}")
                matched.append(("1", fwd, rev))

        matched.sort()

        for i, (rep_id, fwd, rev) in enumerate(matched, start=1):
            key = str(i)
            fwd_dict[key] = fwd
            if not single_end:
                rev_dict[key] = rev

        if not matched:
            raise FileNotFoundError(f"No matching files found for libname '{libname}'")

        lib_entry = {"fwd": fwd_dict}
        if not ignore_umi:
            lib_entry["umi_len"] = umi_len
        if not single_end:
            lib_entry["rev"] = rev_dict
        
        config["libnames"][libname] = lib_entry

    return config

def main():
    opts = argparse.ArgumentParser()

    opts.add_argument("--libnames", dest="libnames", type=list_of_strings, required=True)
    opts.add_argument("--dirs", dest="dirs", type=list_of_strings, required=True)
    opts.add_argument("--single_end", action="store_true", help="Use for single-end reads")
    opts.add_argument("--ignore_umi", action="store_true", help="Use to skip UMI extraction")

    args = opts.parse_args()

    config = make_config(args.libnames, args.dirs, single_end=args.single_end, ignore_umi=args.ignore_umi)

    with open("input-config.json", "w") as outfile:
        json.dump(config, outfile, indent=4)

if __name__ == "__main__":
    main()