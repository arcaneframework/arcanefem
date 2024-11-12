import json
import sys
import os

def find_key(data, target_key):
    """Recursively search for a key in nested dictionaries or lists."""
    if isinstance(data, dict):
        for key, value in data.items():
            if key == target_key:
                return value
            result = find_key(value, target_key)
            if result is not None:
                return result
    elif isinstance(data, list):
        for idx, item in enumerate(data):
            if isinstance(item, dict):
                result = find_key(item, target_key)
                if result is not None:
                    return result
            elif isinstance(item, str):
                if item == target_key and len(data) >= idx + 2:
                    return data[idx + 1]

def load_metrics(config_path):
    """Load metrics from a configuration file."""
    try:
        with open(config_path, 'r') as file:
            return [line.strip() for line in file if line.strip()]
    except IOError as e:
        print(f"Error reading configuration file: {e}")
        sys.exit(1)

def main(file_path, metrics=None, config_path=None):
    formats = [
        "Legacy", "Coo", "CooSort", "Coo_Gpu", "CooSort_Gpu", 
        "Csr", "Csr_Gpu", "CsrNodeWise", "CsrBuildLess"
    ]

    """If there is a format in the metrics given by the user, we only
    search for these specific formats; otherwise, we print for all formats
    present in the JSON file."""
    formats = [fmt for fmt in metrics if fmt in formats] if any(fmt in metrics for fmt in formats) else formats
    
    try:
        with open(file_path, 'r') as file:
            content = file.read()
            obj = json.loads(content)
    except (IOError, json.JSONDecodeError) as e:
        print(f"Error reading or parsing file: {e}")
        sys.exit(1)

    output = []

    output.append('{0: <30}'.format('Number of parallel instance:') + '{0: >5}'.format(str(obj["nbParallelInstance"])))

    cacheWarming = obj["cacheWarming"]

    output.append('{0: <30}'.format('Cache warming:')  + '{0: >5}'.format(str(cacheWarming)))
    
    if 'acceleratorRuntime' in obj:
        output.append('{0: <30}'.format('Accelerator runtime:') + '\"' + '{0: >5}'.format(str(obj["acceleratorRuntime"]) + '\"'))

    output.append('') # newline
    output.append('{0: <30}'.format(f'Dimension:') + '{0: >5}'.format(str(obj["meshDim"])))
    output.append('{0: <30}'.format(f'Node:') + '{0: >5}'.format(str(obj["nbNode"])))
    output.append('{0: <30}'.format(f'Boundary element:') + '{0: >5}'.format(str(obj["nbBoundaryElement"])))
    output.append('{0: <30}'.format(f'Element:') + '{0: >5}'.format(str(obj["nbElement"])) + '\n')

    for format in formats:
        metric_name = f"AssembleBilinearOperator_{format}"
        format_raw_obj = find_key(obj, metric_name)

        if format_raw_obj is not None:
            time = format_raw_obj['Cumulative'].split(' ')[0]

            if cacheWarming > 1:
                time = float(time) / (cacheWarming - 1)

            output.append('{0: <50}'.format(f"{metric_name}:") + f"{time}")

            for key in metrics:
                key_obj = find_key(format_raw_obj, key)

                if key_obj is not None:
                    time = key_obj['Cumulative'].split(' ')[0]

                    if cacheWarming > 1:
                        time = float(time) / (cacheWarming - 1)

                    output.append('{0: <50}'.format(f"  {key}_{format}:") + f"  {time}")
                else:
                    output.append('{0: <50}'.format(f"  {key}_{format}:") + f"  undefined")

            output.append('') # newline between formats

    output.pop() # remove last newline

    for line in output:
        print(line)

if __name__ == "__main__":
    metrics = ["BuildMatrix", "AddAndCompute"] # default sub action metrics to print

    if len(sys.argv) < 2:
        print("Usage: python script.py <json_file_path> <metrics_comma_separated | config_file_path>")
        sys.exit(1)

    file_path = sys.argv[1]

    if len(sys.argv) > 2:
        # Check for command line metrics or config file path
        if len(sys.argv) > 2 and ',' in sys.argv[2] or not os.path.exists(sys.argv[2]):
            metrics = sys.argv[2].split(',')
        else:
            metrics = load_metrics(sys.argv[2])

    main(file_path, metrics)
