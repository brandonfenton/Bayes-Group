## Heatmap with meaningful ordering

# Order the dataset, using the arrange function from the dplyr package
require(dplyr)
od.ordered <- arrange(od, vty, blk, ino, nit)

# Create a matrix of the arranged responses
infected.arranged <- matrix(od.ordered$infected, ncol = 6, byrow = TRUE)

# Set up two panels, right one for a legend
layout(t(1:2), widths = c(9, 1))

# Plot the heatmap, with zeros in black and segments separating the blocks
par(mar = c(3, 10, 6, 2)) # Set big margins
image(z = infected.arranged, y = 1:6,

      # Variety blocks are 1 unit wide, centered at 0.5, 1.5, etc
      x = seq(0.5, 5.5, 1/6),

      # Use black for 0, and use heatmap colors (red-orange-yellow-white) for 1 to 30
      col = c("black", heat.colors(30)), zlim = c(0, 30),

      # Don't automatically create axes or labels
      xlab = "", ylab = "", yaxt = "n", xaxt = "n")

# Use white line segments to visually separate the varieties
segments(x0 = 1.5:4.5, y0 = 0.5, y1 = 6.5, col = "white")

# Place a title at the top, and label Varieties, treatment:status levels, and blocks around the image
title("Infection Counts", line = 4)
axis(3, labels = rep(levels(od.ordered$blk), 5), cex.axis = 0.75,

     # Each block is plotted in a column with width 1/6, so put the labels in the middle
     at = seq(7/12, 5 + 5/12, 1/6))
axis(2, labels = levels(with(od.ordered, interaction(nit, ino))), at = 1:6, las = 2)
axis(1, labels = levels(od.ordered$vty), at = 1:5)

# Legend
par(mar = c(3, 1, 6, 2))
image(y = seq(-0.5, 30.5, 1), z = matrix(0:30, nrow = 1), axes = FALSE, ylab = "",
      col = c("#000000", heat.colors(30)), zlim = c(0, 30))
title("Legend", line = 1.5)
axis(4)
