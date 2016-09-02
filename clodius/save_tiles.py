import collections as col
import cStringIO as csio
import gzip
import itertools as it
import json
import os
import os.path as op
import random
import requests
import slugid
import sys
import time
import itertools

from time import gmtime, strftime

class TileSaver(object):
    def __init__(self, max_data_in_sparse, bins_per_dimension, num_dimensions):
        self.max_data_in_sparse = max_data_in_sparse

        #self.max_data_in_sparse = 0
        self.bins_per_dimension = bins_per_dimension
        self.num_dimensions = num_dimensions

        pass

    def save_tile(self, tile):
        return

    def make_and_save_tile(self, zoom_level, tile_position, tile_data):
        # this implementation shouldn't do anything
        # derived classes should implement this functionality themselves

        #print "saving tile:", zoom_level, tile_position
        tile_id = "{}.{}".format(zoom_level, ".".join(map(str, tile_position)))

        #print "saving:", tile_id
        tile = {'tile_id':  tile_id, "tile_value": tile_data}

        self.save_tile(tile)

    def save_dense_tile(self, zoom_level, tile_position, tile_bins, 
            min_value, max_value):
        initial_values = [0.0] * (self.bins_per_dimension ** self.num_dimensions)

        for (bin_pos, val) in tile_bins.items():
            index = sum([bp * self.bins_per_dimension ** i for i,bp in enumerate(bin_pos)])
            initial_values[index] = val

        self.make_and_save_tile(zoom_level, tile_position, {"dense": 
            [round(v, 5) for v in initial_values],
            'min_value': min_value, 'max_value': max_value })

    def save_sparse_tile(self, zoom_level, tile_position, tile_bins, 
            min_value, max_value):
        shown = []
        for (bin_pos, bin_val) in tile_bins.items():
            shown += [[map(float, bin_pos), bin_val]]

        self.make_and_save_tile(zoom_level, tile_position, {"sparse": shown,
            'min_value': min_value, 'max_value': max_value })

    def save_binned_tile(self, zoom_level, tile_position, tile_bins):
        max_value = max(tile_bins.values())
        min_value = min(tile_bins.values())

        if len(tile_bins) < self.max_data_in_sparse:
            self.save_sparse_tile(zoom_level, tile_position, tile_bins, 
                                  min_value=min_value, max_value=max_value)
        else:
            self.save_dense_tile(zoom_level, tile_position, tile_bins,
                                 min_value=min_value, max_value=max_value)

    def flush():
        return

class EmptyTileSaver(TileSaver):
    def __init__(self, max_data_in_sparse, bins_per_dimension, num_dimensions):
        super(EmptyTileSaver, self).__init__(max_data_in_sparse, 
                                             bins_per_dimension,
                                             num_dimensions)

class ColumnFileTileSaver(TileSaver):
    def __init__(self, max_data_in_sparse, bins_per_dimension, num_dimensions,
            file_path, log_file):
        super(ColumnFileTileSaver, self).__init__(max_data_in_sparse, 
                                             bins_per_dimension,
                                             num_dimensions)
        self.file_path = file_path
        self.bulk_txt = csio.StringIO()
        self.bulk_txt_len = 0
        self.log_file = log_file

        print "created tilesaver:", self.bulk_txt

    def save_tile(self, val):

        '''
        if ('dense' in val['tile_value']):
            value_pos = col.defaultdict(list)
            dense_values = val['tile_value']['dense']
            dense_values = [(x,len(list(y))) for (x,y) in it.groupby(dense_values)]
            dense_values = [item for sublist in dense_values for item in sublist]
            val['tile_value']['dense'] = dense_values
            for i,value in enumerate(dense_values):
                value_pos[value] += [i]
            for key in value_pos:
                sorted_value_pos = sorted(value_pos[key])
                diffs = []
                diffs += [sorted_value_pos[0]]
                for i in range(len(sorted_value_pos)-1):
                    diffs += [sorted_value_pos[i+1] - sorted_value_pos[i]]

                value_pos[key] = diffs
            val['tile_value']['dense'] = value_pos.items()
        '''

        '''
        if ('sparse' in val['tile_value']):
            sparse_values = val['tile_value']['sparse']
            value_pos = col.defaultdict(list)
            for sparse_value in sparse_values:
                value_pos[sparse_value[1]] += [sparse_value[0]]
            val['tile_value']['sparse'] = value_pos.items()

            value_xs_ys = []
            for value, poss in value_pos.items():
                poss = sorted(poss)
                xs = [p[0] for p in poss]
                ys = [p[1] for p in poss]
                value_xs_ys += [value, xs, ys]
            val['tile_value']['sparse'] = value_xs_ys
        '''

        # [[1.0, [[78.0, 123.0], [64.0, 153.0]]]]

        if val["tile_id"] is "tileset_info":
            self.bulk_txt.write(val["tile_id"] + "\t" + "1" + "\t" + "1" + "\t")
        else:
            ti = val['tile_id'].split(".")
            self.bulk_txt.write(str(int(ti[0])+1) + "\t" + str(int(ti[1])+1) + "\t" + str(int(ti[1])+1) + "\t")

        self.bulk_txt.write(json.dumps(val) + "\n")
        curr_pos = self.bulk_txt.tell()
        #print "curr_pos:", curr_pos,self.bulk_txt.getvalue()
        #self.bulk_txt.write(new_string)
        if curr_pos > 2000000:
            self.flush()

    def flush(self):
        if self.bulk_txt.tell() > 0:
            try:
                with open(self.file_path, "a") as column_file:
                    column_file.write(self.bulk_txt.getvalue())
            except Exception as ex:
                if self.log_file is not None:
                    with open(log_file, 'a') as f:
                        f.write(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
                        f.write(ex)

        self.bulk_txt_len = 0
        self.bulk_txt.close()
        self.bulk_txt = csio.StringIO()


class ElasticSearchTileSaver(TileSaver):
    def __init__(self, max_data_in_sparse, bins_per_dimension, num_dimensions,
            es_path, log_file):
        super(ElasticSearchTileSaver, self).__init__(max_data_in_sparse, 
                                             bins_per_dimension,
                                             num_dimensions)
        self.es_path = es_path
        self.bulk_txt = csio.StringIO()
        self.bulk_txt_len = 0
        self.log_file = log_file

        print "created tilesaver:", self.bulk_txt

    def save_tile(self, val):
        # this implementation shouldn't do anything
        # derived classes should implement this functionality themselves

        #self.bulk_txt.write(json.dumps({"index": {"_id": val['tile_id']}}) + "\n")
        if ('sparse' in val['tile_value']):
            sparse_values = val['tile_value']['sparse']

            value_pos = col.defaultdict(list)
            for sparse_value in sparse_values:
                value_pos[sparse_value[1]] += [sparse_value[0]]
            #val['tile_value']['sparse'] = value_pos.items()

            value_xs_ys = []
            for value, poss in value_pos.items():
                poss = sorted(poss)
                dim_values = []
                value_xs_ys += [float(value)]
                value_xs_ys += [float(len(poss))]

                for i in range(len(poss[0])):
                    value_xs_ys += [p[i] for p in poss]

            val['tile_value']['sparse'] = value_xs_ys

        '''
        if ('dense' in val['tile_value']):
            print val['tile_id'], len([x for x in val['tile_value']['dense'] if x > 0])
        '''

        #print "writing:", val['tile_id']
        #val['tile_value']['dense'] = []

        self.bulk_txt.write('{{"index": {{"_id": "{}"}}}}\n'.format(val['tile_id']))
        self.bulk_txt.write(json.dumps(val) + "\n")

        '''
        self.bulk_txt.write('{{"tile_id": {}, "tile_value": '.format(val['tile_id']))

        if 'sparse' in val['tile_value']:
            self.bulk_txt.write(' {{ "sparse": [ ')
            self.bulk_txt.write(','.join(
                ['{{ "pos": [{}], "value": {} }}'.format(
                    ",".join([str(y) for y in x['pos']]), x['value'])
                    for x in val['tile_value']['sparse']]))
            self.bulk_txt.write('] }}')


        #sys.exit(1)
        #new_string += str(val) + "\n"

        #self.bulk_txt_len += len(new_string)
        '''

        curr_pos = self.bulk_txt.tell()
        #self.bulk_txt.write(new_string)
        if curr_pos > 5000000:
            self.flush()

    def flush(self):
        if self.bulk_txt.tell() > 0:
            # only save the tile if it had enough data
            try:
                save_to_elasticsearch("http://" + self.es_path + "/_bulk", self.bulk_txt.getvalue())
            except Exception as ex:
                if self.log_file is not None:
                    with open(log_file, 'a') as f:
                        f.write(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
                        f.write(ex)

            self.bulk_txt_len = 0
            self.bulk_txt.close()
            self.bulk_txt = csio.StringIO()

def save_tile_to_elasticsearch(partition, elasticsearch_nodes, elasticsearch_path):
    bulk_txt = ""
    es_url = op.join(elasticsearch_nodes, elasticsearch_path)
    put_url =  op.join(es_url, "_bulk")

    for val in partition:
        bulk_txt += json.dumps({"index": {"_id": val['tile_id']}}) + "\n"
        bulk_txt += json.dumps(val) + "\n"

        if len(bulk_txt) > 5000000:
            save_to_elasticsearch("http://" + put_url, bulk_txt)
            bulk_txt = ""

    if len(bulk_txt) > 0:
        save_to_elasticsearch("http://" + put_url, bulk_txt)

def save_to_elasticsearch(url, data):
    '''
    Save some data to elastic search.

    The data should be a string suitable for bulk import by
    elasticsearch. The url should be the location of the index, document
    type, along with the _bulk destination.

    :param url: The elasticsearch url that will ingest the data
                e.g. localhost:9200/hg19/tiles/_bulk
    :param data: The data to import.
                e.g. {"index": {"_id": "blah", "my_json": {"x": 2}}}
    '''
    saved = False
    to_sleep = 1
    uid = slugid.nice()
    while not saved:
        try:
            r = requests.post(url, data=data, timeout=8)
            print "Saved", uid,  r, "len(data):", len(data), url #, r.text
            saved = True
            #print "data:", data
        except Exception as ex:

            to_sleep *= 2
            print >>sys.stderr, "Error saving to elastic search (", uid, "), sleeping:", to_sleep, ex
            time.sleep(to_sleep)

            if to_sleep > 600:
                print >>sys.stderr, "Slept too long, returning"
                raise

def save_tile(tile, output_dir, gzip_output):
    '''
    Save a tile to a particular base directory.

    This function create the appropriate sub-directory based on the
    key. 

    They key should be in the format (zoom_level, pos1, pos2...)
    e.g. (5,4,5)
    '''
    key = tile[0]
    tile_value = tile[1]

    outpath = op.join(output_dir, '.'.join(map(str, key)))
    outdir = op.dirname(outpath)

    if not op.exists(outdir):
        try:
            os.makedirs(outdir)
        except OSError as oe:
            # somebody probably made the directory in between when we
            # checked if it exists and when we're making it
            print >>sys.stderr, "Error making directories:", oe

    output_json = {"_source": {"tile_id": ".".join(map(str, key)),
                               "tile_value": tile_value}}
    if gzip_output:
        with gzip.open(outpath + ".gz", 'w') as f:
            f.write(json.dumps(output_json, indent=2))
    else:
        with open(outpath, 'w') as f:
            f.write(json.dumps(output_json, indent=2))
