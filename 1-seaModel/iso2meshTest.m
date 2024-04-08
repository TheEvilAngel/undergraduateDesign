[node1,face1,elem1]=meshabox([0 0 0],[10 10 10],1,1);
[node2,face2,elem2]=meshabox([0 0 0]+5,[10 10 10]+5,1,1);
[newnode,newface]=surfboolean(node1,face1,'union',node2,face2);
plotmesh(newnode,newface);
[newnode,newface]=surfboolean(node1,face1,'-',node2,face2);
figure; plotmesh(newnode,newface);
[newnode,newface]=surfboolean(node1,face1,'diff',node2,face2);
figure;plotmesh(newnode,newface,'x>5');