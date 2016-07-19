c(b=exp_mrna_Mg, c=exp_mrna_Na, a=exp_mrna_Carb)

grid.newpage();
jpeg('rplot.jpg')
example.list = list(A=1:10, B=6:15, C=c(10, 16:20))
venn.grid = venn.diagram(example.list, filename=NULL, print.mode="percent")
grid.draw(venn.grid)
dev.off()


grid.newpage();
exp_mrna_list = list("Carbon \nsource"=exp_mrna_Carb, "Mg stress"=exp_mrna_Mg, "Na stress"=exp_mrna_Na)
venn.grid = venn.diagram(x=exp_mrna_list, filename="exp_mrna_venn.jpeg", 
                         print.mode=c("raw","percent"),force.unique=TRUE, 
                         fill=c("purple","cyan","orange"), 
                         fontface = "bold", cex=rep(1.3,7), cat.fontface="bold", 
                         cat.cex=c(1.4,1.4,1.4),cat.dist=c(0.1, 0.08, 0.05) ,margin = 0.04)
grid.draw(venn.grid)
dev.off()

draw.triple.venn()