from dotenv import load_dotenv
import os


def find_env(var_name: str) -> str:
    load_dotenv()
    value = os.getenv(var_name)
    if value is None:
        raise KeyError(f"Required environment variable '{var_name}' is not set.")
    return value


def find_env_dir(var_name: str) -> str:
    dir_path = find_env(var_name)
    os.makedirs(dir_path, exist_ok=True)

    return dir_path
