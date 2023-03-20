import edlib


def find_all(que, ref, offset=0, max_ed=5, min_length=10):
    results = []
    if len(ref) < min_length:
        return results
    r = edlib.align(que, ref, task="locations", mode="HW")
    ed = r["editDistance"]
    if ed <= max_ed:
        locs = []
        for x, y in r["locations"]:
            y += 1
            if len(locs) == 0:
                locs.append([x, y])
            else:
                if x >= locs[-1][1]:
                    locs.append([x, y])
        for x, y in locs:
            results.append([offset + x, offset + y, ed])
            
        y0 = 0
        data = []
        for x, y in locs:
            if x - y0 >= min_length:
                data.append([ref[y0:x], offset + y0])
            y0 = y
        if len(ref) - y0 >= min_length:
            data.append([ref[y0:len(ref)], offset + y0])
        for item in data:
            for x, y, ed in find_all(que, item[0], item[1], max_ed, min_length):
                results.append([x, y, ed])
    return results