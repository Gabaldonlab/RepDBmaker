import re
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process CD-HIT output clusters")
    parser.add_argument("-i","--input", dest="input_file", required=True, help="Path to the input CD-HIT cluster file")
    parser.add_argument("-o","--output", dest="output_file", required=True, help="Path to the output MMseqs formatted file")
    return parser.parse_args()

def main():
    args = parse_arguments()

    dict_clstr = {}

    with open(args.input_file, "r") as cfh:
        for line in cfh:
            if line.startswith(">Cluster"):
                is_first = True
            elif is_first:
                key = re.findall(r">(.*)\.{3}", line)[0]
                is_first = False
                dict_clstr[key] = [key]
            else:
                other_ids = re.findall(r">(.*)\.{3}", line)[0]
                dict_clstr[key].append(other_ids)

    with open(args.output_file, "w") as cfh:
        for key, values in dict_clstr.items():
            for value in values:
                cfh.write('%s\t%s\n' % (key, value))

if __name__ == "__main__":
    main()
