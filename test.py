from nanonisTCP import nanonisTCP
from nanonisTCP.Bias import Bias

TCP_IP  = '127.0.0.1'                               # Local host
TCP_PORT= 6501                                    # Check available ports in NANONIS > File > Settings Options > TCP Programming Interface
version = 14000                                     # Nanonis RT Engine version number

NTCP = nanonisTCP(TCP_IP, TCP_PORT, version=version)  # This is how you establish a TCP connection. NTCP is the connection handle.
                                                    # Check your nanonis version in Nanonis > help > info and enter the RT Engine number

bias = Bias(NTCP)                                   # Nanonis Bias Module - Pass in the connection handle

bias.Set(1.1)                                       # Set bias to 1.1 V
v = bias.Get()                                      # Get the current bias
print("Bias: " + str(v) + " V")                     # Confirm bias has been set

NTCP.close_connection()                             # Close the connection.

# test