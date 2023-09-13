#!/usr/bin/env python
import optparse


def main():
    
    parser = optparse.OptionParser()
    parser.add_option("-e", "--edit-distance", dest="ed", type="int", default=5, metavar="INT", 
                      help="Maximum of edit distance. [%default]")
    parser.add_option("-p", "--preset", dest="preset", default="SE", metavar="STR", 
                      help="Preset how to split reads. [%default]")
    parser.add_option("-c", "--config", dest="config", metavar="PATH", 
                      help="")
    parser.add_option("-o", "--outdir", dest="outdir", default="./", metavar="STR", 
                      help="")
    options, args = parser.parse_args()
    
    ## TODO
    
    raise NotImplementedError()
    

if __name__ == "__main__":
    main()