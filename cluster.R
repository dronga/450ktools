color.dendrogram.labels = function(n){
      if(is.leaf(n)){
                a = attributes(n)
                        color = foo.lookup.1$color[ which(foo.lookup.1$leaf.name==a$label)[1] ]
                        if( foo.show.labels.1 ){
                                      cex.val = 1
                                                  lab.color = color
                                    }
                        else{
                                      cex.val=0.01
                                                  lab.color="white"
                                    }
                        attr(n, "nodePar") = c(a$nodePar, list(lab.col = lab.color, lab.cex=cex.val, col=color, pch=15, cex=1 ) )
              }
          n
    }

plot.cluster.colors=function( M, colors, force.flat=F, show.labels=T ,main=""){
      ## Hack: sets global environment variable
      d = (dist(t(M)))
          lookup = data.frame( leaf.name = colnames(M), color=colors, stringsAsFactors=F )
          if(force.flat)
                    hangval = -1
          else
                    hangval = 0.1
          D = as.dendrogram(hclust(d), hang=hangval)
          assign("foo.lookup.1", lookup, envir=.GlobalEnv)
          assign("foo.show.labels.1", show.labels, envir=.GlobalEnv)
          dend.colored = dendrapply(D, color.dendrogram.labels)
          plot(dend.colored,main=main)
    }

