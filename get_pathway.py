fin = open('E:\swang141\project\NIH\Data\Stuart Nature ChemBio\Input\\nci_pathway_name.txt')
path = []
path.append([])
for line in fin:
    path.append(line.strip())
fin.close()

fout = open('..\..\embedding\genetic.networktop.pathwayname','w')
fin = open('..\..\embedding\genetic.networktop.pathway')
for line in fin:
    i = int(line.strip())
    fout.write(path[i]+'\n')
fout.close()
