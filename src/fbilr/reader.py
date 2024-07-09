#!/usr/bin/env python

from pygz import PigzFile


class Hit(object):
    def __init__(self, name, direction, location, start, end, ed):
        self.name = name
        self.direction = direction
        self.location = location
        self.start = start
        self.end = end
        self.ed = ed
        
    def __str__(self) -> str:
        return "%s\t%s\t%s\t%d\t%d\t%d" % (self.name, self.direction, self.location, 
                                           self.start, self.end, self.ed)
        

class MatrixRecord(object):
    def __init__(self, name, length, hits, read=None):
        self.name = name
        self.length = length
        self.hits = hits
        self.read = read
    
    def __str__(self):
        items = [self.name, self.length]
        for hit in self.hits:
            items.append(hit.name)
            items.append(hit.direction)
            items.append(hit.location)
            items.append(hit.start)
            items.append(hit.end)
            items.append(hit.ed)
        # if self.read is not None:
        #     pass
        return "\t".join(map(str, items))


class MatrixReader(object):
    def open(path):
        if path.endswith(".gz"):
            f = PigzFile(path)
        else:
            f = open(path)

        for line in f:
            
            row = line.strip("\n").split("\t")
            
            name, length = row[0], int(row[1])
            
            hits = []
            n, m = int((len(row) - 2) / 6), (len(row) - 2) % 6
            assert n > 0
            assert m == 0 or m == 4
            for i in range(n):
                bc, direction, loc, start, end, ed = row[i * 6 + 2:(i + 1) * 6 + 2]
                hit = Hit(name=bc, 
                            direction=direction, 
                            location=loc, 
                            start=int(start), 
                            end=int(end), 
                            ed=int(ed))
                hits.append(hit)
                
            read = None
            if m == 4:
                raise NotImplementedError()
                
            record = MatrixRecord(name=name, length=length, hits=hits, read=read)
                
            yield record
        
        f.close()
        