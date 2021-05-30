# gff3-parser
A sample and powerful gff3 file parser.
# 最好用的GFF3解析工具

## 解析流程

1. 先使用正则表达式提取出第九列中的id。然后存入sqlite。这里有一个小技巧，首先对gff3进行分块，每一行的第三列是其feature的类型，每一个区块都是从`gene`这个类型开始，然后再结束。存储id的区块顺序和id名，这样就可以通过区块来直接进行内容的搜索
2. 给定一个tuple类型的ids，首先根据数据库里面他们对应的区块顺序，进行排序，从而达到对gff文件进行一次遍历，就可以把所有内容全部hit到。
3. 通过`BlockNode`可以实现只对需要的区块解析，不需要的区块可以直接跳过，从而可以大大提高检索的速度。


## example

### 1. 根据id，提取对应的gene所对应的所有feature的行
```python
from gff3 import GFF3
file_path = './test.gff3'
ids = ('rna-XM_008348769.3','gene-LOC114826672')
gff3 = GFF3(file_path,ids)
gff3.search().to_gff3('output.gff3')
```

### 2. 根据id，提取对应的gene所对应的mRNA的行
```python
from gff3 import GFF3, TypeList
file_path = './test.gff3'
ids = ('rna-XM_008348769.3','gene-LOC114826672')
gff3 = GFF3(file_path,ids)
gff3.search(type=TypeList.mrna).to_gff3('output.gff3')
```

### 3. 根据id，同时提取mRNA和cds的注释信息
```python
from gff3 import GFF3, TypeList
file_path = './test.gff3'
ids = ('rna-XM_008348769.3','gene-LOC114826672')
gff3 = GFF3(file_path,ids)
gff3.search(type=TypeList.mrna).to_gff3('output1.gff3')
gff3.search(type=TypeList.cds).to_gff3('output2.gff3')
```

### 4. 使用处理器，将提取出来`start`和`end`转换成整型数字
```python
from gff3 import GFF3, RowNode
file_path = './test.gff3'
ids = ('rna-XM_008348769.3','gene-LOC114826672')

class MyRowNode(RowNode):
    def handle_start(self,start):
        return int(start)
    def handle_end(self,end):
        return int(end)
    
gff3 = GFF3(file_path,ids,rowNodeCls=MyRowNode)
gff3.search().to_json('output.json')
```

### 5. 使用处理器，将提取出来id的前缀去掉
```python
from gff3 import GFF3, RowNode
file_path = './test.gff3'
ids = ('rna-XM_008348769.3','gene-LOC114826672')

class MyRowNode(RowNode):
    def handle_id(self,id):
        if '-' in id:
            return id.split('-')[1]
        return id

gff3 = GFF3(file_path,ids,rowNodeCls=MyRowNode)
gff3.search().to_json('output.json')
```

### 6. 使用全局的处理器，将提取出来的值全部变成大写
```python
from gff3 import GFF3,Handles
file_path = './test.gff3'
ids = ('rna-XM_008348769.3','gene-LOC114826672')
gff3 = GFF3(file_path,ids)

handles = Handles()
handles.add_handle(Handles.upper)  # 自带了一个大写的处理器

gff3.search(handles = handles).to_gff3('output.gff3')
```

[github](https://github.com/DENGWENPENG/gff3-parser)
