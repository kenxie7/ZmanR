zmanr_theme <- theme_pubr() + theme(legend.position="right")+
        theme(
              panel.grid.major.x = element_blank(), #element_line(size = 0.01, colour = "#EAE7E7"),
              panel.grid.major.y = element_blank(), #element_line(size = 0.01, colour = "#EAE7E7"),
              title = element_text(size = 13, color = "black"),
              axis.text.x = element_text(size = 11, color = "black"),
              axis.title.x = element_text(size = 13, color = "black"),#,vjust = -1),
              legend.text = element_text(size = 11, color = "black"),
              legend.title = element_text(size = 13, color = "black"),#,vjust = -1),
              axis.ticks.length.y = unit(.25, "cm"),
              axis.ticks.y = element_line(size = 0.5, color = "black"),
              axis.text.y = element_text(size = 11, color = "black"),
              axis.title.y.left = element_text(size = 13, face = "plain", vjust = 2.2))#,
              #legend.key = element_rect(size = 5, fill="white", colour = NA), legend.key.size = unit(0.50, "cm"), legend.text = element_text(size = 12, face = "plain"),
              #legend.position = c(1.25,0.5) ) #+
              #scale_fill_manual(values =  c("#aa0094","#306674")) + #"#aa0094","#306674"
              #scale_colour_manual(values = c("black","darkgrey")) + #"#aa0094","#306674"
              #scale_y_continuous(limits = c(0, 1), breaks = seq(0,1, by = 0.2)) +
              #scale_x_continuous(limits = c(0, 1), breaks = seq(0,1, by = 0.2)) #+ scale_color_gradient(low="lightgrey", high="red",limits=c(0, 1.0))
