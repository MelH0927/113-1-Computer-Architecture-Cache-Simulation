# 113-1-Computer-Architecture-Cache-Simulation
> 11310CS410000 &lt;Computer Architecture> final project: cache simulation

Run by executing the following lines:

```sh
g++ project.cpp -o project
./project cache.org reference.lst
```

Results will be stored in created file index.rpt.

### Contents of cache.org
```
Address_bits: <num>
Block_size: <num>
Cache_sets: <num>
Associativity: <num
```

### Contents of reference.lst
```
.benchmark <testcase name>
<multiple lines of binary addresses with lenghth = <Address_bits> >
.end
```
