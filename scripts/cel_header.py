import re, struct
from datetime import datetime
def parse_header(payload):
    if len(payload)<24:raise ValueError('Incomplete CEL prefix')
    magic,version=struct.unpack('<ii',payload[:8])
    if (magic,version)!=(64,4):raise ValueError('Unsupported CEL format')
    length=struct.unpack('<i',payload[20:24])[0]
    if not 0<length<100000 or len(payload)<24+length:raise ValueError('Incomplete CEL header')
    header=payload[24:24+length].decode('ascii',errors='strict')
    dat=re.search(r'^DatHeader=(.*)$',header,re.M)
    if not dat:raise ValueError('No DAT header')
    found=re.search(r'\s(\d{2}/\d{2}/\d{2})\s+(\d{2}:\d{2}:\d{2})\s+(\S+)\s+(\S+)',dat.group(1))
    if not found:raise ValueError('Unrecognized acquisition timestamp')
    timestamp=datetime.strptime(found[1]+' '+found[2],'%m/%d/%y %H:%M:%S')
    return dict(array_scan_date=timestamp.date().isoformat(),array_scan_time=timestamp.time().isoformat(),
                scanner_id=found[3],scanner_type=found[4],cel_header=header)

