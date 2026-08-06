# Example usage:
# For DC operating point (OP)
# python3 compare_results.py <custom_op.txt> <ngspice_op.txt> op
# For DC sweep / transient / AC
# python3 compare_results.py <custom_op.txt> <ngspice_op.txt>

import re
import math

def parse_file(filename, type=None):
    data = {}
    numeric_rows = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if type == "op":
                # ngspice
                m = re.match(r'([vViIlL]?\(?[^\s=)]+\)?(?:#branch)?)\s*=\s*([-+0-9.eE]+)', line)
                if m:
                    name, val = m.groups()
                    name = name.lower().replace("#branch", "")
                    match = re.search(r"v\((.*?)\)", name)
                    if match:
                        name = match.group(1)
                    data[name] = float(val)
                    continue

                # custom
                parts = line.strip().split()
                if len(parts) != 2:
                    # skip malformed lines
                    continue
                name, val = parts
                try:
                    data[name.lower()] = float(val)
                except ValueError:
                    # skip lines that don't have numeric values
                    continue
            else:
                try:
                    numeric = [float(x) for x in line.split()]
                    numeric_rows.append(numeric)
                    continue
                except ValueError:
                    pass

    if data:  # OP result
        return data
    elif numeric_rows:  # DC or transient result
        return numeric_rows
    else:
        raise ValueError(f"File '{filename}' contains no recognizable data.")

def mean(values):
    return sum(values) / len(values) if values else 0.0

def compare_results(custom_file, ngspice_file, type=None, tolerance=0.01):
    print(f"Comparing {custom_file} with {ngspice_file}")
    cdata = parse_file(custom_file, type)
    ndata = parse_file(ngspice_file, type)

    if isinstance(cdata, dict) and isinstance(ndata, dict):
        common_keys = sorted(set(cdata.keys()) & set(ndata.keys()))
        if not common_keys:
            print("No common entries found.")
            return

        abs_errs = []
        rel_errs = []
        pass_count = 0

        print(f"{'Name':<15} {'Custom':>12} {'Ngspice':>12} {'Abs Err':>12} {'Rel Err (%)':>12}")
        print("-" * 65)

        for name in common_keys:
            v1 = cdata[name]
            v2 = ndata[name]
            abs_err = abs(v1 - v2)
            rel_err = abs_err / (abs(v2) + 1e-15)
            abs_errs.append(abs_err)
            rel_errs.append(rel_err)
            if rel_err <= tolerance:
                pass_count += 1

            print(f"{name:<15} {v1:12.6e} {v2:12.6e} {abs_err:12.6e} {rel_err*100:12.3f}")

        mae = mean(abs_errs)
        mre = mean(rel_errs)
        pass_rate = pass_count / len(common_keys) * 100

        print("\nSummary:")
        print(f"Mean Absolute Error : {mae:.6e}")
        print(f"Mean Relative Error : {mre*100:.3f}%")
        print(f"Pass Rate (≤{tolerance*100:.1f}%): {pass_rate:.1f}%")

    elif isinstance(cdata, list) and isinstance(ndata, list):
        # Require same shape
        if len(cdata) != len(ndata) or any(len(crow) != len(nrow) for crow, nrow in zip(cdata, ndata)):
            print(f"Shape mismatch: custom {len(cdata)}x{len(cdata[0]) if cdata else 0}, "
                  f"ngspice {len(ndata)}x{len(ndata[0]) if ndata else 0}")
            return

        total_abs = 0.0
        total_rel = 0.0
        total_count = 0
        pass_count = 0

        for crow, nrow in zip(cdata, ndata):
            for v1, v2 in zip(crow, nrow):
                abs_err = abs(v1 - v2)
                rel_err = abs_err / (abs(v2) + 1e-15)
                total_abs += abs_err
                total_rel += rel_err
                total_count += 1
                if rel_err <= tolerance:
                    pass_count += 1

        mae = total_abs / total_count
        mre = total_rel / total_count
        pass_rate = 100 * pass_count / total_count

        print(f"Mean Absolute Error : {mae:.6e}")
        print(f"Mean Relative Error : {mre*100:.3f}%")
        print(f"Pass Rate (≤{tolerance*100:.1f}%): {pass_rate:.1f}%")

    else:
        print("Unsupported data type combination.")
        print(f"custom: {type(cdata)}, ngspice: {type(ndata)}")

if __name__ == "__main__":
    import sys

    custom_file = sys.argv[1]
    ngspice_file = sys.argv[2]
    type = sys.argv[3] if len(sys.argv) > 3 else None

    compare_results(custom_file, ngspice_file, type)
