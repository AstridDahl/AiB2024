# processor name
processor = platform.processor()

# system name
system = platform.system() + " " + platform.release()

# machine architecture
architecture = platform.machine()

# amount of RAM in bytes
ram_bytes = psutil.virtual_memory().total

# Convert bytes to gigabytes (GB)
ram_gb = round(ram_bytes / (1024**3), 2)

# Python version
python_version = platform.python_version()

# Print the system information
print(f"Processor: {processor}")
print(f"System: {system}")
print(f"Architecture: {architecture}")
print(f"RAM: {ram_gb} GB")
print(f"Python version: {python_version}")