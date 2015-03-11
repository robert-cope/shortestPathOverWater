# author: Robert Cope (2015) robert.cope@adelaide.edu.au

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#######################################################################################
#######################################################################################

#This code was developed for the following study:

# Cope RC, Prowse TAA, Ross JV, Wittmann TA, Cassey P (2015) Temporal modelling of ballast water discharge and ship-mediated invasion risk to Australia. Royal Society Open Science. In press.


#######################################################################################
#######################################################################################



#This is an *approximate* shortest path method, over water, between pairs of ports. 
#It will provide a good ballpark in terms of avoiding large land masses, but can get a little distracted around the edges of these land masses.
#if you want a port that isn't in the list, probably the easiest tactic would be to pick a nearby node and use that

library(ncf)
library(geosphere)
library(igraph)
#Csardi G, Nepusz T (2006): The igraph software package for complex network research, InterJournal, Complex Systems 1695. http://igraph.org

#given the nodes and edges, create the graph object
combEdgeList2 <-read.csv('links3combined.txt', header=F,row.names=NULL, stringsAsFactors=F, na.strings='NA',fileEncoding="ISO-8859-1",sep=' ')
nodeList2 <-read.csv('nodes3combined.txt', header=F,row.names=NULL, stringsAsFactors=F, na.strings='NA',fileEncoding="ISO-8859-1",sep=' ')
#the node list contains: (1) grid points in the ocean, (2) nodes at each 'corner' of a continent/island bordering the ocean, and (3) a collection of port locations, by name, from Australia and elsewhere
#edges connect these nodes to nearby nodes
#Antarctica, the top half of siberia, Greenland, and some generally very small islands have been removed to reduce complexity
#based on shorline data from GSHHG http://www.soest.hawaii.edu/pwessel/gshhg/
#these data are described in:
#Wessel, P., and W. H. F. Smith (1996), A Global Self-consistent, Hierarchical, High-resolution Shoreline Database, J. Geophys. Res., 101: 8741-8743.



names(combEdgeList2)[1]<-'node1'
names(combEdgeList2)[2]<-'node2'
names(combEdgeList2)[3]<-'dist'
names(nodeList2)[1]<-'node'
names(nodeList2)[2]<-'long'
names(nodeList2)[3]<-'lat'
nodeList2[,4]<-NULL
combEdgeList2[,4]<-NULL
testGraph2<-graph.edgelist(as.matrix(combEdgeList2[,1:2]),directed=FALSE)
E(testGraph2)$weight=combEdgeList2[,3]
cbind( get.edgelist(testGraph2) , round( E(testGraph2)$weight, 3 ))

#the function that actually gets the paths
getPath<-function(start,end){
  pTest<-get.shortest.paths(testGraph2,which(V(testGraph2)$name==start),which(V(testGraph2)$name==end))
  pN<-V(testGraph2)$name[pTest$vpath[1][[1]]]
  newPN<-pN[1]
  c<-0
  for(i in 2:length(pN)){
    lp<-strsplit(pN[i-1],'')[[1]][1]
    l<-strsplit(pN[i],'')[[1]][1]
    
    if(l=='g' & !lp=='g'){
      #first gridpoint
      newPN<-c(newPN,pN[i])
      c<-0
    } else if(!l=='g' & lp=='g') {
      #last gridpoint
      newPN<-c(newPN,pN[i-1])
      newPN<-c(newPN,pN[i])
    } else if(!l=='g' & !lp=='g'){
      #neither - non-gridpoint
      newPN<-c(newPN,pN[i])
    } else if(l=='g' & lp=='g'){
      y2<-nodeList2[nodeList2$node==pN[i],][[3]]
      ln<-strsplit(pN[i+1],'')[[1]][1]
      if(ln=='g' & y2 > -50){
        adjList<-V(testGraph2)[nei(which(V(testGraph2)$name==pN[i]))]
        for(k in 1:length(adjList$name)){
          nodeA<-adjList$name[k]
          fl<-strsplit(nodeA,'')
          if (fl[[1]][1]=='w'){
            newPN<-c(newPN,nodeA)
            break
          }
        }
      }
    }
  }
  newPN2<-data.frame(node=character(),long=numeric(),
                     lat=numeric(),
                     stringsAsFactors=FALSE)
  newPN2<-rbind(newPN2,c(newPN[1],nodeList2[nodeList2$node == newPN[1],2:3]))
  names(newPN2)[1]<-'node'
  newPN2 <- transform(newPN2, node = as.character(node))
  for(i in 2:length(newPN)){
    lp<-strsplit(newPN[i-1],'')[[1]][1]
    l<-strsplit(newPN[i],'')[[1]][1]
    if(!l=='g') {
      newPN2<-rbind(newPN2,c(newPN[i],nodeList2[nodeList2$node == newPN[i],][2][[1]],nodeList2[nodeList2$node == newPN[i],][3][[1]]))
    } else if(l=='g' & lp=='g'){
      p1<-nodeList2[nodeList2$node == newPN[i-1],2:3]
      p2<-nodeList2[nodeList2$node == newPN[i],2:3]
      dGC<-gcdist(p1[1][[1]],p1[2][[1]],p2[1][[1]],p2[2][[1]])
      path<-gcIntermediate(p1,p2,ceiling(dGC/100))
      for (j in 1:length(path[,1])){
        newPN2<-rbind(newPN2,c(as.character('int'),path[j,][[1]],path[j,][[2]]))
      }
    }
  }
  newPN2$pol<-rep(-1,length(newPN2$long))
  newPN2$long<-as.numeric(newPN2$long)
  newPN2$lat<-as.numeric(newPN2$lat)
  return(newPN2)
}


getPathLengthGC<-function(path){
  l<-0
  for (i in 1:(length(path$long)-1)){
    l<-l+gcdist(path[i,]$long,path[i,]$lat,path[i+1,]$long,path[i+1,]$lat)
  }
  return(l)
}

#e.g.
msTestPath<-getPath('PORT-MORESBY','TOWNSVILLE')
msTestPath<-getPath('FREMANTLE','ROTTERDAM')
msTestPath<-getPath('MELBOURNE','SINGAPORE')
msTestPath<-getPath('DAMPIER','QINGDAO')
msTestPath<-getPath('NEWCASTLE','QINGDAO')
msTestPath<-getPath('FREMANTLE','RIO-DE-JANEIRO-RJ')
msTestPath<-getPath('BRISBANE','GUAM')
msTestPath<-getPath('FREMANTLE','ROTTERDAM')

getPathLengthGC(msTestPath)

