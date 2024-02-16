

g_bind <- ggarrange(g1, g2 + labs(y = "", fill = "Variation", color = "Variation"), 
          g3 + theme(legend.position = "none"), 
          g4 + labs(y = "")  + theme(legend.position = "none"), 
          align = "v", heights = c(1,0.8), labels = "AUTO", label.x = 0.1, label.y = c(0.85,0.85, 1, 1))
ggsave("graphs/figure2.png",g_bind, dpi = 600, width = 6.5, height = 6.5)