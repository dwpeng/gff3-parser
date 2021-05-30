# coding:utf-8
# @author:邓文鹏
# @email:1732889554@qq.com
# @update_time:2021/5/29
# @desc:This code segment is based this document that from https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md


import re
import sqlite3
import os
import json

# ID匹配pattern
ID_DEFAULT_PATTERN = '\tID=([\w\.\-]+);?'


def set_id_pattern(patt):
    global ID_DEFAULT_PATTERN
    ID_DEFAULT_PATTERN = patt


def reset_id_pattern():
    global ID_DEFAULT_PATTERN
    ID_DEFAULT_PATTERN = '\tID=([\w\.\-]+);?'


class GffNodeType(str):
    ...


class TypeList:
    # The types of the feature
    gene = GffNodeType('gene')
    mrna = GffNodeType('mrna')
    cds = GffNodeType('cds')
    lnc_rna = GffNodeType('lnc_rna')
    exon = GffNodeType('exon')
    region = GffNodeType('region')
    pseudogene = GffNodeType('pseudogene')
    transcript = GffNodeType('transcript')
    unknown = GffNodeType('unknown')


class Handles:

    def __init__(self, handles: tuple = None):
        if handles and isinstance(handles, tuple):
            self.__handles = list(handles)
        self.__handles = []

    @staticmethod
    def lower(x):
        return str(x).lower()

    @staticmethod
    def upper(x):
        return str(x).upper()

    def add_handle(self, handle):
        self.__handles.append(handle)
        return self

    def pop_handle(self):
        self.__handles.pop()
        return self

    def handles(self):
        return self.__handles

    def __str__(self):
        return f"[{', '.join(['Handle(' + i.__name__  + ',' + str(index + 1) +')' for index,i in enumerate(self.__handles)])}]"

    def do_handles(self, raw):
        result = [raw] + [handle(raw) for handle in self.__handles]
        return result[-1]


class Handle_result:
    def __init__(self, raw, handle):
        self.handle = handle
        self.raw = raw
        self.result = None

    def __call__(self):
        if self.handle:
            self.result = self.handle(self.raw)
        self.result = self.raw
        return self.result


class RowNode:
    def __init__(self, raw, handles: Handles):
        self._raw = raw
        if handles:
            if not isinstance(handles, Handles):
                raise TypeError(
                    'A Handles type is needed, but got a %s type ' % type(handles))
            self.raw = handles.do_handles(raw)  # 使用处理器处理
        else:
            self.raw = raw
        self.split_raw = None
        vaild_result = self._vaild_gff3_row(self.raw)
        if vaild_result:
            self.split_raw = vaild_result
        else:
            raise ValueError()
        # 根据每一列的位置，进行解析
        self.seqid = self.handle_seqid(self.split_raw[0])
        self.source = self.handle_source(self.split_raw[1])
        self.type = self.handle_type(self.split_raw[2])
        self.start = self.handle_start(self.split_raw[3])
        self.end = self.handle_end(self.split_raw[4])
        self.score = self.handle_score(self.split_raw[5])
        self.strand = self.handle_strand(self.split_raw[6])
        self.phase = self.handle_phase(self.split_raw[7])
        # 属性单独的解析逻辑
        self.attrs = self.handle_attrs(self._parser_attrs, self.split_raw[8])
        # 对gene id和parent属性进行单独处理，方便后面进行搜索的操作
        self.id = self.handle_id(self.attrs.get('id', None))
        self.parent = self.handle_parent(self.attrs.get('parent', None))

        # self.is_gene = 'gene' in str(self.type).lower()
        # self.gene_uid = uuid.uuid4()  # 如果是gene，就生成一个单独的uid，方便对每个基因下的外显子，转录本进行归类

    def _vaild_gff3_row(self, raw: str, sep='\t'):
        split = raw.split(sep)
        result = []
        # 检查一下是否有9行
        split_columns = len(split)
        result.append(split_columns == 9)
        # TODO
        # 这里还可以添加更多的判断条件
        # 后期可以改成可以自动添加判断逻辑的程序
        return split if all(result) else False

    # 完全可以重写的处理器
    def handle_seqid(self, seqid): return seqid
    def handle_source(self, source): return source
    def handle_type(self, type): return type
    def handle_start(self, start): return start
    def handle_end(self, end): return end
    def handle_score(self, score): return score
    def handle_strand(self, strand): return strand
    def handle_phase(self, phase): return phase
    def handle_attrs(self, parser_handle, attrs): return parser_handle(attrs)

    # 处理属性内的id和paerent的值，可以做字典的替换
    def handle_id(self, id):
        return id

    def handle_parent(self, parent):
        return parent

    def _parser_attrs(self, attrs):
        attrs_map = {}
        if ';' in attrs:
            attrs = attrs.split(';')
        else:
            attrs = [attrs]
        for one_attr in attrs:
            k, v = one_attr.split('=')
            attrs_map[k.lower()] = v
        # 最后一定要return出来解析好的attr_map
        return attrs_map

    def to_dict(self):
        return dict(
            raw=self._raw,
            seqid=self.seqid,
            source=self.source,
            type=self.type,
            score=self.score,
            start=self.start,
            end=self.end,
            strand=self.strand,
            phase=self.phase,
            attrs=self.attrs,
            id=self.id,
            parent=self.parent
        )

    @classmethod
    def is_row_node(cls):
        """用来做继承时的类型判断"""
        return True

class BlockNode:
    """这里根据gff3注释的习惯，来进行分块
        type 从gene开始，到gene结束
    """

    def __init__(self, raw):
        self.raw = raw
        self.row_nodes = None

    def to_row_nodes(self, handles: Handles, rowNodeCls=None):
        # 先把每一行进行划分
        raws = self.raw.split('\n')
        # 这里可以实现对rowNode定制功能，比如修改处理器的功能，只需要一个继承的类就可以了
        if not rowNodeCls:
            rowNodeCls = RowNode
        else:
            try:
                is_row_node = rowNodeCls.is_row_node()
                if not is_row_node:
                    raise TypeError()
            except:
                raise TypeError('A RowNode type is expected, but got a %s type.' % type(rowNodeCls))
        self.row_nodes = [rowNodeCls(i, handles)
                          for i in raws if i and '#' not in i]  # 把注释行和空行去掉

    def to_dict(self):
        if not self.row_nodes:
            return []
        return [i.to_dict() for i in self.row_nodes]

    @property
    def keys(self):
        return [i.lower() for i in re.findall(ID_DEFAULT_PATTERN, self.raw, re.S | re.I)]


class DB:
    conn = None

    def create_db(self, db_path, blockNodeList):
        if not os.path.exists(db_path):
            conn = sqlite3.connect(db_path)
            conn.execute('''CREATE TABLE id_index 
            (id INT  NOT NULL,
            id_name TEXT NOT NULL);
            ''')
            c = conn.cursor()
            block_count = 0
            insert = 'INSERT INTO id_index(id,id_name) VALUES '
            for blockNode in blockNodeList:
                sql = ''
                for key in blockNode.keys:
                    sql += '('
                    sql += f'{block_count},"{key}"'  # 区块顺序 和 id名
                    sql += '),'
                block_count += 1
                sql = insert + sql.strip(',')
                c.execute(sql)
            conn.commit()
            self.conn = conn
        else:
            self.conn = sqlite3.connect(db_path)

    def search_and_sort_ids(self, ids):
        """根据建好的数据库，可以查询出各个id的先后顺序，然后进行排序，这样把gff3文件从头过一遍就可以把所有的id全部查出来"""
        if not self.conn:
            raise ValueError('You should run `create_db` first.')
        _ids = []
        for id in ids:
            search_sql = 'select id,id_name from id_index where id_name like "%s"' % str(
                id).lower()
            rows = self.conn.execute(search_sql)
            rows = list(rows)
            if rows:
                _ids.append(rows[0])
        _ids.sort(key=lambda x: x[0], reverse=True)  # 根据区块顺序进行排序，区块越大越靠前
        return _ids


# 初始化一个数据库实例
db = DB()


class GFF3:

    def __init__(self, gff3_file_path, ids: tuple, rowNodeCls=None):
        """
            Parameters:
                - gff3_file_path: gff3文件地址
                - ids: 一个存放ids的tuple
                - rowNodeCls: 来自定义每一行处理的行为，受json中的cls参数启发
        """
        if not os.path.exists(gff3_file_path):
            raise FileExistsError('File not exist.Please check again!')
        self.gff3_file = open(gff3_file_path, 'r')
        db.create_db(gff3_file_path + '.db', self.to_block())
        self.ids = self.read_ids(ids)
        self.rowNodeCls = rowNodeCls

    def to_block(self):
        gff = self.gff3_file
        node = ''
        for line in gff:
            if not line or '#' in line:
                continue
            if f'{TypeList.gene}\t' in line and node:
                yield BlockNode(node)
                node = ''
            node = node + line
        yield BlockNode(node)
        self.gff3_file.seek(0)  # 把文件的指针移动到0，后面还需要从头过一遍来做真正的内容提取

    def read_ids(self, ids):
        return db.search_and_sort_ids(ids)

    def search(self, *, type=None, handles=None):
        block_count = 0
        result = []
        if not self.ids:
            return
        ids = [(0, 0)] + [i for i in self.ids]  # 做一份拷贝
        count, _ = ids.pop()
        for block in self.to_block():
            if count == block_count:
                result.append(block)
                self._check_same_block(block_count, ids)  # 处理相同区块的情况
                count, _ = ids.pop()
            block_count += 1
        if handles and not isinstance(handles, Handles):
            raise TypeError(
                'handles expected a Handle type, but got a %s.' % type(handles))
        [i.to_row_nodes(handles,rowNodeCls=self.rowNodeCls) for i in result]
        result = [i.to_dict() for i in result]
        # 在这里实现type的过滤功能
        if type:
            _result = []
            for block in result:
                _result.append(
                    list(
                        filter(lambda x: x['type'].lower() == type.lower(), block))
                )
            result = _result
        return Result(result)

    def _check_same_block(self, block_count, ids):
        if ids[-1][0] == block_count:
            ids.pop()
            self._check_same_block(block_count, ids)  # 通过递归的方式，确保每个区块都只查找一次


class Result:
    cols = ['seqid', 'source', 'type', 'start',
            'end', 'score', 'strand', 'phase', 'attrs']

    def __init__(self, result):
        self.result = result

    def _save(self, txt, path):
        with open(path, 'w') as f:
            f.write(txt)

    def _to_row(self, line):
        """"把gff的一行从字典转化成文本形式"""
        cols = [line.get(i, '') for i in self.cols]
        if isinstance(cols[-1], dict):
            cols[-1] = ';'.join([k + '=' + v for k, v in cols[-1].items()])
        return '\t'.join(cols)

    def to_gff3(self, path):
        txt = ''
        for r in self.result:
            for line in r:
                txt += self._to_row(line)
                txt += '\n'
            txt += '#\n'
        txt = txt.strip('#\n')
        self._save(txt, path)

    def to_json(self, path):
        txt = json.dumps(self.result)
        self._save(txt, path)

    def __str__(self):
        from pprint import pprint
        pprint(self.result)
        return ''


if __name__ == '__main__':
    ids = ('gene-LOC103431662', 'cds-XP_028955284.1')
    path = 'GFF3/genomic.gff3'
    gff3 = GFF3(path, ids)
    handles = Handles()
    print(gff3.search(type=TypeList.mrna, handles=handles).to_gff3('mrna.gff3'))
    print(gff3.search(type=TypeList.gene, handles=handles).to_gff3('gene.gff3'))
