import json
import sys

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
                    return data[idx + 1];

metrics_to_extract = ["BuildMatrix", "AddAndCompute"]

def extract_metrics(data, keys):
    """Extract 'Local' and 'Cumulative' metrics for specified keys from nested data."""
    extracted_data = {}
    for key in keys:
        key_obj = find_key(data, key)
        if key_obj is not None:
            extracted_data[key] = {
                "Local": key_obj.get("Local"),
                "Cumulative": key_obj.get("Cumulative")
            }
    return extracted_data

def main(file_path):
    formats = [
        "Legacy", "Coo", "CooSort", "Coo_Gpu", "CooSort_Gpu", 
        "Csr", "Csr_Gpu", "CsrNodeWise", "CsrBuildLess"
    ]
    
    try:
        with open(file_path, 'r') as file:
            content = file.read()
            obj = json.loads(content)
    except (IOError, json.JSONDecodeError) as e:
        print(f"Error reading or parsing file: {e}")
        sys.exit(1)

    output = {}
    for format in formats:
        metric_name = f"AssembleBilinearOperator_{format}"
        format_raw_obj = find_key(obj, metric_name)

        if format_raw_obj is not None:
            format_obj = {"Local": format_raw_obj["Local"], "Cumulative": format_raw_obj["Cumulative"]}

            # Extract additional metrics using helper function
            additional_metrics = extract_metrics(format_raw_obj, metrics_to_extract)
            format_obj.update(additional_metrics)

            output[format] = {"AssembleBilinearOperator": format_obj}
            format_raw_obj = find_key(obj, metric_name)

    print(json.dumps(output, indent=4))

if __name__ == "__main__":
    if len(sys.argv) > 1:
        file_path = sys.argv[1]
        main(file_path)
    else:
        print("Usage: python script.py <json_file_path>")

