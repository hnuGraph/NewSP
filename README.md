# NewSP
## Introduction
In this study, we address the problem of unnecessary computations in traditional continuous subgraph matching (CSM) frameworks due to premature expansions of the search space in dynamic graphs. Traditional CSM frameworks expand small partial matches according to a specific matching order until the final results are obtained. This extension involves two sequential steps: computing candidate vertices for an unmapped query vertex and expanding the search space using these candidate data. However, this long-established search model has a potential flaw, as premature expansions of the search space can lead to unnecessary computations. To address this issue, we introduce a novel search process, NewSP. Unlike traditional methods, NewSP emphasizes operations rather than extensions, incorporating a unique feature of postponing expansion at the operation level. This approach prevents premature expansions without compromising the initial pruning power of the selected matching order. Furthermore, NewSP allows for multiple consecutive expansions, paving the way for a multi-expansion strategy for further optimization. Our model also enables the implementation of cache strategies for candidate set reuse, as it does not necessitate immediate expansion of a candidate set once identified. To improve performance, we propose an adaptive index filtering strategy independent of the specific index used. Comprehensive experiments demonstrate that our method improves by up to two to three orders of magnitude compared to traditional algorithms.

## Compile

Our framework requires c++17 and GCC 7.x (or later). One can compile the code by executing the following commands. 

```shell
make
```

## Execute

After a successful compilation, the binary file is created under the `build/` directory. One can execute CSM using the following command.

```shell
build/csm -q <query-graph-path> -d <data-graph-path> -u <update-stream-path> 
```


### Commandline Parameters

Other commandline parameters supported by the framework are listed in the following table.

| Command Line Parameters | Description                                                     | Valid Value      | Default Value |
|-------------------------|-----------------------------------------------------------------|------------------|---------------|
| --time-limit            | Time limit for the incremental matching phase (in seconds).     | 0-4294967295     | 3600          |
| --report-initial        | Perform initial matching or not.                                | on/off           | on            |
| --initial-time-limit    | Time limit for the initial matching phase (in seconds).         | 0-4294967295     | 4294967295    |
| --print-prep            | Print preprocessing results or not.                             | on/off           | on            |
| --print-enum            | Print matches results or not.                                   | on/off           | off           |


For example, if one requires the framework (1) to return after finding the first result on each update operation; and (2) to spend at most 1 hour (3600 seconds) on the incremental matching, then the command should be

```shell
build/csm -q <query-graph-path> -d <data-graph-path> -u <update-stream-path> -a <algorithm> --time-limit 3600
```

## Input File Format
Both the input query graph and data graph are vertex- and edge-labeled. Each vertex is represented by a distinct unsigned integer (from 0 to 4294967295). There is at most one edge between two arbitrary vertices. 

### Query Graph

Each line in the query graph file represent a vertex or an edge.

1. A vertex is represented by `v <vertex-id> <vertex-label>`;
2. An edge is represented by `e <vertex-id-1> <vertex-id-2> <edge-label>`.

The two endpoints of an edge must appear before the edge. For example, 

```
v 0 0
v 1 0
e 0 1 0
v 2 1
e 0 2 1
e 2 1 2
```

### Initial Data Graph

The initial data graph file has the same format as the query graph file.

### Graph Update Stream

Graph update stream is a collection of insertions and deletions of a vertex or an edge.

1. A vertex insertion is represented by `v <vertex-id> <vertex-label>`;
2. A vertex deletion is represented by `-v <vertex-id> <vertex-label>`;
3. An edge insertion is represented by `e <vertex-id-1> <vertex-id-2> <edge-label>`;
4. An edge deletion is represented by `-e <vertex-id-1> <vertex-id-2> <edge-label>`;

The vertex or edge to be deleted must exist in the graph, and the label must be the same as that in the graph. If an edge is inserted to the data graph, both its endpoints must exist. For example,

```
v 3 1
e 2 3 2
-v 2 1
-e 0 1 0
```

##  Datasets

We provide 4 datasets in our experiment 

1. Amazon dataset        [[link](https://snap.stanford.edu/data/com-Amazon.html)]
2. Livejournal dataset   [[link](https://snap.stanford.edu/data/soc-LiveJournal1.html)]
3. LSBench dataset       [[link](https://code.google.com/archive/p/lsbench/)]
4. Netflow dataset       [[link](https://catalog.caida.org/dataset/passive\_2013\_pcap)]



**Summary of Datasets**

| **Datasets** |         **Type**        | **Vertexes** |  **Edges**  | **Average Degree** |
| :----------: |    :-----------------:  | :----------: | :---------: | :----------------: |
|    Amazon    |     Product network     |    403,394   |  2,433,408  |       12.06        |
|  Livejournal |    Community network    |   4,847,571  | 42,841,237  |       17.68        |
|   LSBench    | Benchmark data generator|   5,210,099  | 20,270,676  |       7.78         |
|   Netflow    |     Network traffic     |   3,114,895  |  2,849,732  |       1.83         |
