## kartoon plot for AS visualization
# for reverse arrow direction
# replace scale_x_continuous with scale_x_reverse

SE <- function(reverse=F) {
  if (reverse==F) {
    color1 <- "#4DAF4A"
    color2 <- "#984EA3"
  } else {
    color2 <- "#4DAF4A"
    color1 <- "#984EA3"
  }
  #
  # prepare kartoon data
  # SE, Cassettes
  AS_1 <- data.frame(EXONSTART=c(0,2,4),
                     EXONEND=c(1,3,5),
                     EXONSTRAND="+")
  seg_1 <- data.frame(x=c(0, 1.5, 3.5),
                      xend=c(1.5, 3.5, 4.5),
                      y=0, yend=0)
  AS_2 <- data.frame(EXONSTART=c(0,4),
                     EXONEND=c(1,5),
                     EXONSTRAND="+")
  seg_2 <- data.frame(x=c(0, 1.5, 2.5, 3.5),
                      xend=c(1.5, 2.5, 3.5, 4.5),
                      y=0, yend=0)
  #
  options(repr.plot.width=2.5, repr.plot.height=1)
  p1 <- ggplot() +
    # add line and arrow
    # geom_hline(yintercept=0)+
    # geom_segment(x = 0, xend = 5.5, y = 0, yend = 0) +
    geom_segment(data=seg_1, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length=unit(0.3,"cm")), size=1) +
    # add exons
    geom_rect(data=AS_1, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.1,ymax=0.1),fill=color1) +
    # add UTR
    # geom_rect(data=UTR, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.05,ymax=0.05),fill="#282a73")+
    # details
    labs(title = NULL,subtitle = NULL)+
    theme_void() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.1, 0.1), limits = c(0, 5)) +
    scale_y_continuous(expand = c(0.02, 0.02), limits = c(-0.1, 0.1))
  p1
  #
  options(repr.plot.width=2.5, repr.plot.height=1)
  p2 <- ggplot() +
    # add line and arrow
    # geom_hline(yintercept=0)+
    # geom_segment(x = 0, xend = 5.5, y = 0, yend = 0) +
    geom_segment(data=seg_2, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length=unit(0.3,"cm")), size=1) +
    # add exons
    geom_rect(data=AS_2, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.1,ymax=0.1),fill=color2)+
    # add UTR
    # geom_rect(data=UTR, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.05,ymax=0.05),fill="#282a73")+
    # details
    labs(title = NULL,subtitle = NULL)+
    theme_void() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.1, 0.1), limits = c(0, 5))+
    scale_y_continuous(expand = c(0.02, 0.02), limits = c(-0.1, 0.1))
  p2
  #
  if (reverse==F) {
    options(repr.plot.width=2.5, repr.plot.height=2)
    SE1 <- cowplot::plot_grid(p1,p2,ncol = 1)
    SE1
  } else {
    options(repr.plot.width=2.5, repr.plot.height=2)
    SE2 <- cowplot::plot_grid(p2,p1,ncol = 1)
    SE2
  }
}

AF <- function(reverse=F) {
  #
  if (reverse==F) {
    color1 <- "#4DAF4A"
    color2 <- "#984EA3"
  } else {
    color2 <- "#4DAF4A"
    color1 <- "#984EA3"
  }
  # prepare kartoon data
  # AF, Alternative first exon
  AS_1 <- data.frame(EXONSTART=c(2,4),
                     EXONEND=c(3,5),
                     EXONSTRAND="+")
  seg_1 <- data.frame(x=c(3,3.5),
                      xend=c(3.5,5),
                      y=0, yend=0)
  #
  AS_2 <- data.frame(EXONSTART=c(1,4),
                     EXONEND=c(2,5),
                     EXONSTRAND="+")
  seg_2 <- data.frame(x=c(2,3),
                      xend=c(3,5),
                      y=0, yend=0)
  #
  UTR_1 <- data.frame(EXONSTART=c(1),
                      EXONEND=c(2),
                      EXONSTRAND="+")
  UTR_2 <- data.frame(EXONSTART=c(0),
                      EXONEND=c(1),
                      EXONSTRAND="+")
  #
  options(repr.plot.width=2.5, repr.plot.height=1)
  p1 <- ggplot() +
    # add line and arrow
    # geom_hline(yintercept=0)+
    # geom_segment(x = 0, xend = 5.5, y = 0, yend = 0) +
    geom_segment(data=seg_1, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length=unit(0.3,"cm")), size=1) +
    # add exons
    geom_rect(data=AS_1, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.1,ymax=0.1),fill=color1)+
    # add UTR
    geom_rect(data=UTR_1, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.05,ymax=0.05),fill=color1)+
    # details
    labs(title = NULL,subtitle = NULL)+
    theme_void() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.1, 0.1), limits = c(0, 5)) +
    scale_y_continuous(expand = c(0.02, 0.02), limits = c(-0.1, 0.1))
  p1
  #
  options(repr.plot.width=2.5, repr.plot.height=1)
  p2 <- ggplot() +
    # add line and arrow
    # geom_hline(yintercept=0)+
    # geom_segment(x = 0, xend = 5.5, y = 0, yend = 0) +
    geom_segment(data=seg_2, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length=unit(0.3,"cm")), size=1) +
    # add exons
    geom_rect(data=AS_2, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.1,ymax=0.1),fill=color2)+
    # add UTR
    geom_rect(data=UTR_2, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.05,ymax=0.05),fill=color2)+
    # details
    labs(title = NULL,subtitle = NULL)+
    theme_void() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.1, 0.1), limits = c(0, 5))+
    scale_y_continuous(expand = c(0.02, 0.02), limits = c(-0.1, 0.1))
  p2
  #
  if (reverse==F) {
    options(repr.plot.width=2.5, repr.plot.height=2)
    AF1 <- cowplot::plot_grid(p1,p2,ncol = 1)
    AF1
  } else {
    options(repr.plot.width=2.5, repr.plot.height=2)
    AF2 <- cowplot::plot_grid(p2,p1,ncol = 1)
    AF2
  }
}

AL <- function(reverse = F) {
  #
  if (reverse==F) {
    color1 <- "#4DAF4A"
    color2 <- "#984EA3"
  } else {
    color2 <- "#4DAF4A"
    color1 <- "#984EA3"
  }
  # prepare kartoon data
  # AF, Alternative first exon
  AS_1 <- data.frame(EXONSTART=c(0,2),
                     EXONEND=c(1,3),
                     EXONSTRAND="+")
  seg_1 <- data.frame(x=c(1,1.5),
                      xend=c(1.5,2.5),
                      y=0, yend=0)
  #
  AS_2 <- data.frame(EXONSTART=c(0,3),
                     EXONEND=c(1,4),
                     EXONSTRAND="+")
  seg_2 <- data.frame(x=c(1,2),
                      xend=c(2,3.5),
                      y=0, yend=0)
  #
  UTR_1 <- data.frame(EXONSTART=c(3),
                      EXONEND=c(4),
                      EXONSTRAND="+")
  UTR_2 <- data.frame(EXONSTART=c(4),
                      EXONEND=c(5),
                      EXONSTRAND="+")
  #
  options(repr.plot.width=2.5, repr.plot.height=1)
  p1 <- ggplot() +
    # add line and arrow
    # geom_hline(yintercept=0)+
    # geom_segment(x = 0, xend = 5.5, y = 0, yend = 0) +
    geom_segment(data=seg_1, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length=unit(0.3,"cm")), size=1) +
    # add exons
    geom_rect(data=AS_1, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.1,ymax=0.1),fill=color1)+
    # add UTR
    geom_rect(data=UTR_1, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.05,ymax=0.05),fill=color1)+
    # details
    labs(title = NULL,subtitle = NULL)+
    theme_void() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.1, 0.1), limits = c(0, 5)) +
    scale_y_continuous(expand = c(0.02, 0.02), limits = c(-0.1, 0.1))
  p1
  #
  options(repr.plot.width=2.5, repr.plot.height=1)
  p2 <- ggplot() +
    # add line and arrow
    # geom_hline(yintercept=0)+
    # geom_segment(x = 0, xend = 5.5, y = 0, yend = 0) +
    geom_segment(data=seg_2, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length=unit(0.3,"cm")), size=1) +
    # add exons
    geom_rect(data=AS_2, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.1,ymax=0.1),fill=color2)+
    # add UTR
    geom_rect(data=UTR_2, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.05,ymax=0.05),fill=color2)+
    # details
    labs(title = NULL,subtitle = NULL)+
    theme_void() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.1, 0.1), limits = c(0, 5))+
    scale_y_continuous(expand = c(0.02, 0.02), limits = c(-0.1, 0.1))
  p2
  #
  #
  if (reverse==F) {
    options(repr.plot.width=2.5, repr.plot.height=2)
    AL1 <- cowplot::plot_grid(p1,p2,ncol = 1)
    AL1
  } else {
    options(repr.plot.width=2.5, repr.plot.height=2)
    AL2 <- cowplot::plot_grid(p2,p1,ncol = 1)
    AL2
  }
}

A3 <- function(reverse = F) {
  #
  if (reverse==F) {
    color1 <- "#4DAF4A"
    color2 <- "#984EA3"
  } else {
    color2 <- "#4DAF4A"
    color1 <- "#984EA3"
  }
  # prepare kartoon data
  # A3, Alternative acceptor
  AS_1 <- data.frame(EXONSTART=c(0,3),
                     EXONEND=c(1,5),
                     EXONSTRAND="+")
  seg_1 <- data.frame(x=c(1, 2),
                      xend=c(2, 5),
                      y=0, yend=0)
  AS_2 <- data.frame(EXONSTART=c(0,4),
                     EXONEND=c(1,5),
                     EXONSTRAND="+")
  seg_2 <- data.frame(x=c(1, 2, 3),
                      xend=c(2, 3, 5),
                      y=0, yend=0)
  #
  options(repr.plot.width=2.5, repr.plot.height=1)
  p1 <- ggplot() +
    # add line and arrow
    # geom_hline(yintercept=0)+
    # geom_segment(x = 0, xend = 5.5, y = 0, yend = 0) +
    geom_segment(data=seg_1, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length=unit(0.3,"cm")), size=1) +
    # add exons
    geom_rect(data=AS_1, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.1,ymax=0.1),fill=color1)+
    # add UTR
    # geom_rect(data=UTR, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.05,ymax=0.05),fill="#282a73")+
    # details
    labs(title = NULL,subtitle = NULL)+
    theme_void() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.1, 0.1), limits = c(0, 5)) +
    scale_y_continuous(expand = c(0.02, 0.02), limits = c(-0.1, 0.1))
  p1
  #
  options(repr.plot.width=2.5, repr.plot.height=1)
  p2 <- ggplot() +
    # add line and arrow
    # geom_hline(yintercept=0)+
    # geom_segment(x = 0, xend = 5.5, y = 0, yend = 0) +
    geom_segment(data=seg_2, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length=unit(0.3,"cm")), size=1) +
    # add exons
    geom_rect(data=AS_2, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.1,ymax=0.1),fill=color2)+
    # add UTR
    # geom_rect(data=UTR, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.05,ymax=0.05),fill="#282a73")+
    # details
    labs(title = NULL,subtitle = NULL)+
    theme_void() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.1, 0.1), limits = c(0, 5))+
    scale_y_continuous(expand = c(0.02, 0.02), limits = c(-0.1, 0.1))
  p2
  #
  if (reverse==F) {
    options(repr.plot.width=2.5, repr.plot.height=2)
    A3_1 <- cowplot::plot_grid(p1,p2,ncol = 1)
    A3_1
  } else {
    options(repr.plot.width=2.5, repr.plot.height=2)
    A3_2 <- cowplot::plot_grid(p2,p1,ncol = 1)
    A3_2
  }
}

A5 <- function(reverse = F) {
  #
  if (reverse==F) {
    color1 <- "#4DAF4A"
    color2 <- "#984EA3"
  } else {
    color2 <- "#4DAF4A"
    color1 <- "#984EA3"
  }
  # prepare kartoon data
  # A3, Alternative acceptor
  AS_1 <- data.frame(EXONSTART=c(0,4),
                     EXONEND=c(2,5),
                     EXONSTRAND="+")
  seg_1 <- data.frame(x=c(2, 3),
                      xend=c(3, 5),
                      y=0, yend=0)
  AS_2 <- data.frame(EXONSTART=c(0,4),
                     EXONEND=c(1,5),
                     EXONSTRAND="+")
  seg_2 <- data.frame(x=c(1, 2, 3),
                      xend=c(2, 3, 5),
                      y=0, yend=0)
  #
  options(repr.plot.width=2.5, repr.plot.height=1)
  p1 <- ggplot() +
    # add line and arrow
    # geom_hline(yintercept=0)+
    # geom_segment(x = 0, xend = 5.5, y = 0, yend = 0) +
    geom_segment(data=seg_1, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length=unit(0.3,"cm")), size=1) +
    # add exons
    geom_rect(data=AS_1, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.1,ymax=0.1),fill=color1)+
    # add UTR
    # geom_rect(data=UTR, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.05,ymax=0.05),fill="#282a73")+
    # details
    labs(title = NULL,subtitle = NULL)+
    theme_void() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.1, 0.1), limits = c(0, 5)) +
    scale_y_continuous(expand = c(0.02, 0.02), limits = c(-0.1, 0.1))
  p1
  #
  options(repr.plot.width=2.5, repr.plot.height=1)
  p2 <- ggplot() +
    # add line and arrow
    # geom_hline(yintercept=0)+
    # geom_segment(x = 0, xend = 5.5, y = 0, yend = 0) +
    geom_segment(data=seg_2, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length=unit(0.3,"cm")), size=1) +
    # add exons
    geom_rect(data=AS_2, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.1,ymax=0.1),fill=color2)+
    # add UTR
    # geom_rect(data=UTR, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.05,ymax=0.05),fill="#282a73")+
    # details
    labs(title = NULL,subtitle = NULL)+
    theme_void() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.1, 0.1), limits = c(0, 5))+
    scale_y_continuous(expand = c(0.02, 0.02), limits = c(-0.1, 0.1))
  p2
  #
  if (reverse==F) {
    options(repr.plot.width=2.5, repr.plot.height=2)
    A5_1 <- cowplot::plot_grid(p1,p2,ncol = 1)
    A5_1
  } else {
    options(repr.plot.width=2.5, repr.plot.height=2)
    A5_2 <- cowplot::plot_grid(p2,p1,ncol = 1)
    A5_2
  }
}

MX <- function(reverse = F) {
  #
  if (reverse==F) {
    color1 <- "#4DAF4A"
    color2 <- "#984EA3"
  } else {
    color2 <- "#4DAF4A"
    color1 <- "#984EA3"
  }
  # prepare kartoon data
  # A3, Alternative acceptor
  AS_1 <- data.frame(EXONSTART=c(0, 1.75, 4),
                     EXONEND=c(1, 2.25, 5),
                     EXONSTRAND="+")
  seg_1 <- data.frame(x=c(1, 1.5, 2, 3),
                      xend=c(1.5, 2, 3, 4.5),
                      y=0, yend=0)
  AS_2 <- data.frame(EXONSTART=c(0, 2.75, 4),
                     EXONEND=c(1, 3.25, 5),
                     EXONSTRAND="+")
  seg_2 <- data.frame(x=c(1, 2, 3, 3.75),
                      xend=c(2, 3, 3.75, 4.5),
                      y=0, yend=0)
  #
  options(repr.plot.width=2.5, repr.plot.height=1)
  p1 <- ggplot() +
    # add line and arrow
    # geom_hline(yintercept=0)+
    # geom_segment(x = 0, xend = 5.5, y = 0, yend = 0) +
    geom_segment(data=seg_1, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length=unit(0.3,"cm")), size=1) +
    # add exons
    geom_rect(data=AS_1, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.1,ymax=0.1),fill=color1)+
    # add UTR
    # geom_rect(data=UTR, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.05,ymax=0.05),fill="#282a73")+
    # details
    labs(title = NULL,subtitle = NULL)+
    theme_void() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.1, 0.1), limits = c(0, 5)) +
    # scale_x_reverse() +
    scale_y_continuous(expand = c(0.02, 0.02), limits = c(-0.1, 0.1))
  p1
  #
  options(repr.plot.width=2.5, repr.plot.height=1)
  p2 <- ggplot() +
    # add line and arrow
    # geom_hline(yintercept=0)+
    # geom_segment(x = 0, xend = 5.5, y = 0, yend = 0) +
    geom_segment(data=seg_2, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length=unit(0.3,"cm")), size=1) +
    # add exons
    geom_rect(data=AS_2, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.1,ymax=0.1),fill=color2)+
    # add UTR
    # geom_rect(data=UTR, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.05,ymax=0.05),fill="#282a73")+
    # details
    labs(title = NULL,subtitle = NULL)+
    theme_void() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.1, 0.1), limits = c(0, 5))+
    # scale_x_reverse() +
    scale_y_continuous(expand = c(0.02, 0.02), limits = c(-0.1, 0.1))
  p2
  #
  if (reverse==F) {
    options(repr.plot.width=2.5, repr.plot.height=2)
    MX1 <- cowplot::plot_grid(p1,p2,ncol = 1)
    MX1
  } else {
    options(repr.plot.width=2.5, repr.plot.height=2)
    MX2 <- cowplot::plot_grid(p2,p1,ncol = 1)
    MX2
  }
}

RI <- function(reverse = F) {
  #
  if (reverse==F) {
    color1 <- "#4DAF4A"
    color2 <- "#984EA3"
  } else {
    color2 <- "#4DAF4A"
    color1 <- "#984EA3"
  }
  # prepare kartoon data
  # A3, Alternative acceptor
  AS_1 <- data.frame(EXONSTART=c(0, 4),
                     EXONEND=c(1, 5),
                     EXONSTRAND="+")
  seg_1 <- data.frame(x=c(1, 1.5, 2.5, 3.5, 4.5),
                      xend=c(1.5, 2.5, 3.5, 4.5, 4.75),
                      y=0, yend=0)
  AS_2 <- data.frame(EXONSTART=c(0),
                     EXONEND=c(5),
                     EXONSTRAND="+")
  seg_2 <- data.frame(x=c(1, 2, 3, 3.75),
                      xend=c(2, 3, 3.75, 4.5),
                      y=0, yend=0)
  #
  options(repr.plot.width=2.5, repr.plot.height=1)
  p1 <- ggplot() +
    # add line and arrow
    # geom_hline(yintercept=0)+
    # geom_segment(x = 0, xend = 5.5, y = 0, yend = 0) +
    geom_segment(data=seg_1, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length=unit(0.3,"cm")), size=1) +
    # add exons
    geom_rect(data=AS_1, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.1,ymax=0.1),fill=color1)+
    # add UTR
    # geom_rect(data=UTR, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.05,ymax=0.05),fill="#282a73")+
    # details
    labs(title = NULL,subtitle = NULL)+
    theme_void() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.1, 0.1), limits = c(0, 5)) +
    scale_y_continuous(expand = c(0.02, 0.02), limits = c(-0.1, 0.1))
  p1
  #
  options(repr.plot.width=2.5, repr.plot.height=1)
  p2 <- ggplot() +
    # add line and arrow
    # geom_hline(yintercept=0)+
    # geom_segment(x = 0, xend = 5.5, y = 0, yend = 0) +
    geom_segment(data=seg_2, aes(x=x,y=y,xend=xend,yend=yend), arrow=arrow(length=unit(0.3,"cm")), size=1) +
    # add exons
    geom_rect(data=AS_2, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.1,ymax=0.1),fill=color2)+
    # add UTR
    # geom_rect(data=UTR, aes(xmin=EXONSTART, xmax=EXONEND,ymin= -0.05,ymax=0.05),fill="#282a73")+
    # details
    labs(title = NULL,subtitle = NULL)+
    theme_void() +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0.1, 0.1), limits = c(0, 5))+
    scale_y_continuous(expand = c(0.02, 0.02), limits = c(-0.1, 0.1))
  p2
  #
  if (reverse==F) {
    options(repr.plot.width=2.5, repr.plot.height=2)
    RI1 <- cowplot::plot_grid(p1,p2,ncol = 1)
    RI1
  } else {
    options(repr.plot.width=2.5, repr.plot.height=2)
    RI2 <- cowplot::plot_grid(p2,p1,ncol = 1)
    RI2
  }
}

