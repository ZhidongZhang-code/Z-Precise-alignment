import os
import tkinter as tk
from tkinter import filedialog, messagebox

# 导入您的其他函数和模块

def run_script():
    """
    运行主要的脚本功能
    """
    try:
        main()
        messagebox.showinfo("Success", "Analysis completed successfully!")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {str(e)}")

def browse_file(entry):
    """
    打开文件对话框以选择文件，并将选择的路径显示在文本框中
    """
    filename = filedialog.askopenfilename()
    entry.delete(0, tk.END)
    entry.insert(0, filename)

def browse_folder(entry):
    """
    打开文件夹对话框以选择文件夹，并将选择的路径显示在文本框中
    """
    foldername = filedialog.askdirectory()
    entry.delete(0, tk.END)
    entry.insert(0, foldername)

def create_gui():
    """
    创建图形用户界面
    """
    # 创建主窗口
    root = tk.Tk()
    root.title("Bioinformatics Analysis Tool")

    # 创建输入文件路径选择部件
    input_label = tk.Label(root, text="Input Fasta File:")
    input_label.grid(row=0, column=0, padx=10, pady=5)
    input_entry = tk.Entry(root, width=50)
    input_entry.grid(row=0, column=1, padx=10, pady=5)
    input_button = tk.Button(root, text="Browse", command=lambda: browse_file(input_entry))
    input_button.grid(row=0, column=2, padx=5, pady=5)

    # 创建输出文件路径选择部件
    output_label = tk.Label(root, text="Output Folder:")
    output_label.grid(row=1, column=0, padx=10, pady=5)
    output_entry = tk.Entry(root, width=50)
    output_entry.grid(row=1, column=1, padx=10, pady=5)
    output_button = tk.Button(root, text="Browse", command=lambda: browse_folder(output_entry))
    output_button.grid(row=1, column=2, padx=5, pady=5)

    # 创建运行按钮
    run_button = tk.Button(root, text="Run Analysis", command=run_script)
    run_button.grid(row=2, column=1, padx=10, pady=10)

    # 运行主循环
    root.mainloop()

if __name__ == "__main__":
    create_gui()
