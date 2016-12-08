# Args:
#    network: file name of input
#    Example: karate.txt
#    input network should be an edge list
#############################################
BNC<-function(network){
    library('igraph')
    library('rlist')
   #############################################
   #input the network file
   E.lists <-read.table("/users/yevisliu/Desktop/network",header = FALSE)
   aj.mat<-as.matrix(E.lists)
   label <- paste(aj.mat[, 1], aj.mat[, 2], sep="-")
   # transfer into edge lists of graph
   g<-graph.edgelist(aj.mat,directed=FALSE)
   # return the number of nodes and edges
   n<-vcount(g)
   nodes<-seq(n)
   m<-ecount(g)
   # Truncation threshold
   threshold<-round(sqrt(n))
   # degree of nodes
   d<-degree(g)
   # adjacent list of each node
   r<-get.adjlist(g)
   #############################################
   # k-shell decompose coefficient of nodes
   K.core<-graph.coreness(g)
   # extract 1-core nodes
   n.core1<-which(K.core==1)
   #############################################

   #############################################
   #   take all the roads of this network      #
   #############################################
   # return the shortest paths of nodes in the network % grouped by nodes %
   fun<-function(x){get.all.shortest.paths(g,x)$res}
   sp.lists<-sapply(seq(n),fun)
   fun1<-function(x){sapply(get.all.shortest.paths(g,x)$res,length)}
   sp.length<-sapply(seq(n),fun1)
   sp.list<-list.ungroup(sp.lists)
   sp.le<-unlist(sp.length)
   # count the number of branchs (weights) of the shortest paths from one vertex to the other one % grouped by nodes %
   fun2<-function(x){get.all.shortest.paths(g,x)$nrgeo}
   sp.branch<-sapply(seq(n),fun2)
   sp.br<-split(sp.branch, rep(1:ncol(sp.branch), each = nrow(sp.branch)))
   # add the corresponding weights to each path
   fun3<-function(x){rep(x,x)}
   n.sp.br<-list.ungroup(sapply(sp.br,fun3))
   # return all the paths which length larger than two
   l.path<-sp.list[sp.le>2]
   # return the number of paths between any node pairs
   n.path<-n.sp.br[sp.le>2]
   # return correponding weights of above paths
   w.path<-1/n.path
   #############################################

   #############################################
   # add penalty factor for those roads contain 1-core nodes
   #############################################
   # transfer l.path into matrix
   ls<-lapply(l.path, 'length<-', max(lengths(l.path)))
   lmat<-do.call(rbind,ls)
   lmat[is.na(lmat)]<-0
   fa<-function(x)(is.element(n.core1,x))
   le<-apply(lmat,1,fa)
   ifelse(length(n.core1)==1,le<-t(t(le)),le<-t(le))
   fb<-function(x)(which(x==T))
   ga<-apply(le,2,fb)
   gb<-sort(unique(unlist(ga)))
   # extract corresponding weights
   wa<-w.path[gb]
   wb<-wa/max(d)
   W.sp<-replace(w.path,gb,wb)
   #############################################

   #############################################
   # delete two vertices of each shortest path
   #############################################
   fun4<-function(x){x[-1]}
   fun5<-function(x){x[1:length(x)-1]}
   Pnodes<-sapply(l.path,fun4)
   Pnodes<-sapply(Pnodes,fun5)
   # reorder the elements for each path
   Spnodes<-sapply(Pnodes,sort)
   #############################################

   #############################################
   #      count the weights for each road      #
   #############################################
   n.list<-lapply(Spnodes, 'length<-', max(lengths(Spnodes)))
   p.mat<-do.call(rbind,n.list)
   p.mat[is.na(p.mat)]<-0
   # add weight for each path
   plist<-cbind(p.mat,W.sp)
   colnames(plist)<-c(seq(dim(plist)[2]))
   colnames(plist)[dim(plist)[2]]<-"weight"
   plist<-as.data.frame(plist)
   # count the frequency based on the first four columns
   p.count<-aggregate(weight ~ ., data = plist, sum)
   p.count$weight<-p.count$weight/2
   #############################################
   
   #############################################
   #            Rorder all the roads           #
   #############################################
   nc<-ncol(p.count)
   p.count<-p.count[order(p.count[,nc],decreasing=T),]
   W.sp<-p.count$weight
   p.road<-p.count[-c(nc)]
   p.road[p.road==0]<-NA
   v<-as.list(as.data.frame(t(p.road)))
   u<-lapply(v, function(x) x[!is.na(x)])
   u.length<-sapply(u,length)
   # return all the roads
   u.road<-u[u.length>=2]
   W.road<-W.sp[u.length>=2]
   R.length<-u.length[u.length>=2]
   u.roads<-lapply(u.road,'length<-',max(lengths(u.road)))
   R.mat<-do.call(rbind,u.roads)
   R.mat[is.na(R.mat)]<-0
   # check whether each node is the element of R.mat
   R.node<-as.list(data.frame(t(seq(n))))
   fun6<-function(x){is.element(R.node,x)}
   R.element<-apply(R.mat,1,fun6)
   fun7<-function(x){which(x==T)}
   R.list<-apply(R.element,1,fun7)
   R.lists<-sapply(R.list,list.ungroup)
   # divide the frequency of each road by its length of path
   R.fre<-W.road/(R.length)
   # return the road betweenness of nodes
   fun8<-function(x){sum(R.fre[x])}
   Road.betweenness.node<-sapply(R.lists,fun8)
   #############################################

   #############################################
   #   compute the bridging node centrality    #
   #############################################
   # get all the shortest paths for each pairs of nodes
   sp<-shortest.paths(g)
   sp[is.infinite(sp)]<-0
   sp.sum<-rowSums(sp)
   ##################################
   # count the number of those paths which length=1 or 2
   fc<-function(x)(length(x[x==1]))
   fd<-function(x)(length(x[x==2]))
   num1<-apply(sp,1,fc)
   num2<-apply(sp,1,fd)
   sp.sum<-sp.sum-(num1+num2*2)
   ##################################
   # road betweeness
   bnc.betweenness<-Road.betweenness.node/(d)
   bnc.betweenness[is.na(bnc.betweenness)]<-0
   bnc.betweenness<-(bnc.betweenness-min(bnc.betweenness))/(max(bnc.betweenness)-min(bnc.betweenness))
   # The bridging coefficient
   bnc.coefficient<-1/(sp.sum)
   bnc.coefficient[is.infinite(bnc.coefficient)]<-0
   bnc.coefficient<-(bnc.coefficient-min(bnc.coefficient))/(max(bnc.coefficient)-min(bnc.coefficient))
   # The bridging centrality
   bnc.centrality<-bnc.coefficient*bnc.betweenness
   names(bnc.centrality)<-seq(n)
   Me<-cbind(bnc.centrality,nodes)
   index<-Me[order(Me[, 1], decreasing=T), ]
   rownames(index)<-seq(n)
   # print the top results
   print(head(index, threshold))
}
#############################################
