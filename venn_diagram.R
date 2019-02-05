library(venneuler)
v <- venneuler(c(Depth=0, Distance=1, BottomTemp=0, BottomSal=2,
                 "Depth&BottomTemp"= 1))
plot(v)

library(VennDiagram)
library(wesanderson)

col.palette <- wes_palette("Darjeeling2", 4, type = "discrete")
palette(col.palette)

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local Adaptation/venn_diag_blank.png", width=5, height=5, res=300, units="in")

par(
  mar=c(5, 4, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14 # point size, which is the font size
)

venn.plot <- draw.quad.venn(
  area1 = 9,
  area2 = 6,
  area3 = 6,
  area4 = 10,
  n12 = 2,
  n13 = 5,
  n14 = 9,
  n23 = 2,
  n24 = 2,
  n34 = 5,
  n123 = 1,
  n124 = 2,
  n134 = 5,
  n234 = 1,
  n1234 = 1,
  # category = c("Distance", "Bottom \nSalinity", "Depth", "Bottom \nTemperature" ),
  # fill = c("orange", "tomato", "green", "blue"),
  # category =c c("", "", "", "" ),
  fill = col.palette,
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("black", "black", "black", "black"),
  cat.dist = c(.21,.24,.1,.12)
)

dev.off()

