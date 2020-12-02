import gzip
import os
from to_10x import ToTenx
from tqdm import tqdm
import argparse
import  sys
from functools import partial

def parse_args():
    description = "input data directory to be transform and output dicrectory!"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input_dir', '-i', type=str, default='./totenxtest/datadir', help='the input data directory!')
    parser.add_argument('--output_dir', '-o', type=str, default='./totenxtest/outdatadir', help='the output data directory!')
    parser.add_argument('--trim', '-trm', type=str, default='"', help='tirm the symbol! quo for \'  and dque for " ')
    parser.add_argument('--seperater', '-spr', type=str, default=' ', help='tirm the symbol! s for space ')

    args = parser.parse_args()
    return args

def batchtransform(source_dir, out_dir, datapro=None):
    temfile = 'temporarytenxfile.txt'
    if os.path.exists(temfile):
        os.remove(temfile)
    allfiles = os.listdir(source_dir)
    print('{num} files was found!'.format(num=str(len(allfiles))))
    for eachfile in tqdm(allfiles):
        if os.path.isfile(os.path.join(source_dir, eachfile)):
            filename, _type = os.path.splitext(eachfile)
            out_dir_each = os.path.join(out_dir, filename)
            try:
                os.mkdir(out_dir_each)
            except Exception as e:
                print(e)
                print("directory alread exists!")
                sys.exit()
            if os.path.exists(temfile):
                os.remove(temfile)

            if datapro:
                if _type =='.gz':
                    content = gzip.GzipFile(os.path.join(source_dir, eachfile))
                    content = content.read().decode('utf-8')
                    content = datapro(content)
                    open(temfile, 'wb+').write(content.encode('utf-8'))
                elif _type =='.txt':
                    content = open(os.path.join(source_dir, eachfile), 'r', encoding='utf-8').read()
                    content = datapro(content)
                    open(temfile, 'wb+').write(content.encode('utf-8'))
                else:
                    print("error! unrecognize file format!")
                    sys.exit()

                totenx = ToTenx(temfile, out_dir_each)
            else:
                if _type =='.gz':
                    content = gzip.GzipFile(os.path.join(source_dir, eachfile))
                    open(temfile, 'wb+').write(content.read())
                    totenx = ToTenx(temfile, out_dir_each)
                else:
                    totenx = ToTenx(os.path.join(source_dir, eachfile), out_dir_each)


    if os.path.exists(temfile):
        os.remove(temfile)
def data_process(content, trim_sig='', sepr= ''):
    if trim_sig:
        content = content.replace(trim_sig, '')
    if sepr:
        content = content.replace(sepr, '\t')
    return content



if __name__ == '__main__':
    args = parse_args()
    if args.trim == 'que':
        args.trim = "'"
    if args.trim == 'dque':
        args.trim = '"'
    if args.seperater == 's':
        args.seperater = ' '


    if args.trim or args.seperater:
        data_process = partial(data_process, trim_sig=args.trim, sepr=args.seperater)
        batchtransform(args.input_dir, args.output_dir, datapro=data_process)
    else:
        batchtransform(args.input_dir, args.output_dir)
    #input_dir = './totenxtest/datadir'
    # output_dir = './totenxtest/outdatadir'
    # batchtransform(input_dir, output_dir)
