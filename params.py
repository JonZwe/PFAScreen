import yaml
import os
import pandas as pd

def get_config(config_path="parameters.yaml"):
    """
    Load parameters from YAML file and extract sample names
    
    Args:
        config_path: Path to YAML config file
        
    Returns:
        Tuple of (config_dict, sample_names, output_folder)
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        params = yaml.safe_load(f)
    
    # Extract sample names from sample file
    sample_df = pd.read_csv(params['path_sample_file'])
    sample_names = sample_df['sample'].tolist()
    
    # Extract output folder
    output_folder = params['path_output']
    
    return params, sample_names, output_folder


if __name__ == "__main__":
    # Test the function
    try:
        config, sample_names, output_folder = get_config()
        print("Loaded parameters:")
        for key, value in config.items():
            print(f"  {key}: {value}")
        print(f"\nSample names: {sample_names}")
        print(f"Output folder: {output_folder}")
    except Exception as e:
        print(f"Error: {e}")
