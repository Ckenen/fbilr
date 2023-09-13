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
                pass
                
            record = MatrixRecord(name=name, length=length, hits=hits, read=read)
                
            yield record
        
        f.close()
        

class Matrix2Record(object):
    def __init__(self):
        # read
        self.read_name = None
        self.read_length = None
        # head
        self.head_barcode = "."
        self.head_direction = "."
        self.head_location = "."
        self.head_start = -1
        self.head_end = -1
        self.head_ed = -1
        # tail
        self.tail_barcode = "."
        self.tail_direction = "."
        self.tail_location = "."
        self.tail_start = -1
        self.tail_end = -1
        self.tail_ed = -1
        
    def __str__(self):
        row = [self.read_name, self.read_length, 
               self.head_barcode, self.head_direction, self.head_location, 
               self.head_start, self.head_end, self.head_ed, 
               self.tail_barcode, self.tail_direction, self.tail_location, 
               self.tail_start, self.tail_end, self.tail_ed]
        s = "\t".join(map(str, row))
        return s


class Matrix2Reader(object):
    @classmethod
    def open(cls, path):
        with PigzFile(path, "rt") as f:
            for line in f:
                row = line.strip("\n").split("\t")
                record = Matrix2Record()
                record.read_name = row[0]
                record.read_length = int(row[1])
                if row[2] != ".":
                    record.head_barcode = row[2]
                    record.head_direction = row[3]
                    record.head_location = row[4]
                    record.head_start = int(row[5])
                    record.head_end = int(row[6])
                    record.head_ed = int(row[7])
                if row[8] != ".":
                    record.tail_barcode = row[8]
                    record.tail_direction = row[9]
                    record.tail_location = row[10]
                    record.tail_start = int(row[11])
                    record.tail_end = int(row[12])
                    record.tail_ed = int(row[13])
                yield record
