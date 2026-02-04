import socket

host = "127.0.0.1"
port = 6501

with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
    s.settimeout(2)  # Timeout in Sekunden
    try:
        s.connect((host, port))
        print(f"✅ Port {port} ist offen")
    except (ConnectionRefusedError, socket.timeout):
        print(f"❌ Port {port} ist geschlossen")
    # close the connection automatically with 'with' statement

    #s.close()
    #print(f"Verbindung zu {host}:{port} geschlossen")