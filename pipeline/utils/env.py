from dotenv import load_dotenv
import os


def find_env(var_name: str) -> str:
    load_dotenv()
    value = os.getenv(var_name)
    if value is None:
        raise KeyError(f"Required environment variable '{var_name}' is not set.")
    return value


def find_env_dir(var_name: str) -> str:
    root_path = find_env("ROOT_DIR")
    dir_path = find_env(var_name)
    full_path = os.path.join(root_path, dir_path)
    os.makedirs(full_path, exist_ok=True)

    return full_path
