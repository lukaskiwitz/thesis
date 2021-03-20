import os
import glob
from shutil import copyfile, rmtree

from parameters import path


for collect, format in [("slice","png"),("volume","png")]:
    imgdir = os.path.join(path, "images")

    scans = {}

    for scan_dir in glob.glob1(imgdir,"scan_*"):
        scan_index = scan_dir

        scan_dir = os.path.join(imgdir,scan_dir)
        fields = {}
        for field_dir in glob.glob1(scan_dir,"*"):
            field_name = field_dir
            if field_name not in fields.keys():
                fields[field_name] = []
            field_dir = os.path.join(scan_dir, field_dir)
            for ts_dir in glob.glob1(field_dir,"*"):
                ti = ts_dir
                ts_dir = os.path.join(field_dir,ts_dir)
                img_path = os.path.join(ts_dir,collect+"."+format)
                try:
                    fields[field_name].append((int(ti),img_path))
                except:
                    print("")
        scans[scan_index] = fields


    imgdir = os.path.join(imgdir,collect)
    if os.path.isdir(imgdir):
        rmtree(imgdir)

    for scan_index, fields in scans.items():
        scan_img_dir = os.path.join(imgdir,scan_index)
        for field_name,file_list in fields.items():
            field_dir = os.path.join(scan_img_dir,field_name)
            os.makedirs(field_dir,exist_ok=True)
            for ti, file_path in file_list:
                copyfile(file_path,os.path.join(field_dir,"{f}_{ti}.{form}".format(ti = ti,f = field_name,form = format)))


