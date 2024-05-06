#set margin parametres for the plot
par(mar=c(3,3,3,3))

#make the empty plot
plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 10), ylim = c(0, 10),
     xaxs = "i", yaxs = "i", axes = F, frame.plot = TRUE)

#make lines creating a table (cols)
lines(c(0, 0), c(0, 10), lwd = 3)
lines(c(1, 1), c(0, 10), lwd = 2)
lines(c(5.5, 5.5), c(0, 10))
lines(c(10, 10), c(0, 10), lwd = 3)

#make lines creating a table (rows)
lines(c(0, 10), c(0, 0), lwd = 3)
lines(c(1, 10), c(10/7, 10/7))
lines(c(1, 10), c(10/7 *2, 10/7 *2))
lines(c(0, 10), c(10/7 *3, 10/7 *3), lwd = 2)
lines(c(1, 10), c(10/7 *4, 10/7 *4))
lines(c(1, 10), c(10/7 *5, 10/7 *5))
lines(c(0, 10), c(10/7 *6, 10/7 *6), lwd = 2)
lines(c(0, 10), c(10, 10), lwd = 3)

#add text to the table
text(0.5, 2.2, 'Precipitation', srt = 90, cex = 1.5, font = 2)
text(0.5, 6.3, 'Temperature', srt = 90, cex = 1.5, font = 2)

text(3.2, 0.7, 'max PPT', cex = 1.5)
text(3.2, 10/7 + 0.7, 'mean PPT', cex = 1.5)
text(3.2, (10/7 * 2) + 0.7, 'max PPT', cex = 1.5)
text(3.2, (10/7 * 3) + 0.7, 'min T', cex = 1.5)
text(3.2, (10/7 * 4) + 0.7, 'mean T', cex = 1.5)
text(3.2, (10/7 * 5) + 0.7, 'max T', cex = 1.5)

text(7.7, 0.7, 'mean T', cex = 1.5)
text(7.7, 10/7 + 0.7, 'mean T', cex = 1.5)
text(7.7, (10/7 * 2) + 0.7, 'mean T', cex = 1.5)
text(7.7, (10/7 * 3) + 0.7, 'mean PPT', cex = 1.5)
text(7.7, (10/7 * 4) + 0.7, 'mean PPT', cex = 1.5)
text(7.7, (10/7 * 5) + 0.7, 'mean PPT', cex = 1.5)

text(3.2, (10/7 * 6) + 0.7, 'Dependend variable', cex = 1.5, font = 2)
text(7.7, (10/7 * 6) + 0.7, 'Control variable', cex = 1.5, font = 2)



