import datetime

def print_message(*message):
    timestamp = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    print("[" + timestamp + "] " + " ".join(map(str, message)))
