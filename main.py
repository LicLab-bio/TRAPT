import argparse
import json
import os
import sys

if __name__=='__main__':
    try:
        parser = argparse.ArgumentParser(description='')
        parser.add_argument('-i','--input',dest='input',required=True,type=str,help='Gene list of the targets. Only official gene symbol.')
        parser.add_argument('-o','--output',dest='output',required=True,type=str,help='Output file path.')
        parser.add_argument('--iter',dest='vgae_iterations',type=int,default=10, required=False,help='Vgae iteration times.')
        parser.add_argument('--threshold',dest='range_threshold',type=float,default=0.05, required=False,help='High specificity enhancer signal threshold.')
        args = parser.parse_args()
        if not os.path.exists(args.output):
            os.mkdir(args.output)
        with open(os.path.join(args.output,"params.json"),"w") as file:
            file.write(json.dumps(args.__dict__))
        os.system(f"bash run.sh {args.output}")
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted me!\n")
        sys.exit(0)
