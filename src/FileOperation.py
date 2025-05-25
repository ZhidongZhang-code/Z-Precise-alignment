import os

def remove_file(file):
    """
    删除指定文件。

    Args:
        file (str): 要删除的文件路径。
    """
    if os.path.exists(file):
        os.remove(file)

# 加载YAML配置文件的函数
def load_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)