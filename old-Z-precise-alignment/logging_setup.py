import argparse
from logging_setup import setup_logging

def main():
    parser = argparse.ArgumentParser(description='Bioinformatics tool script.')
    parser.add_argument('--log', type=str, help='Log file path (optional)')
    parser.add_argument('--log_level', type=str, default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help='Set the log level (default: INFO)')
    args = parser.parse_args()

    # 设置日志
    logger = setup_logging(log_file_path=args.log, log_level=args.log_level.upper())

    # 使用 logger 记录消息
    logger.info("Script started")
    # 这里进行其他操作...
    logger.info("Script finished")

if __name__ == "__main__":
    main()
