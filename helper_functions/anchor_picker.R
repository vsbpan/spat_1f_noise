library(shiny)
library(tidyverse)

anchor_picker_app <- function(img_path, 
                            thin = 3,
                            anchor_size = 3){
  
  ui_anchor_picker <- fluidPage(
    shiny::titlePanel(
      htmlOutput("scoring_status")
    ),
    wellPanel(
      actionButton("remove_anchors", "Reset"),
      actionButton("reproject_button", "Reproject"),
      actionButton("save", "Save score")
    ),
    shiny::mainPanel(
      imageOutput("plot", dblclick = "image_dblclick")
    )
  )
  
  
  
  server_anchor_picker <- function(input, output, session) {
    
    output$scoring_status <- renderText({
      HTML(paste0("<b>", gsub(".*/","",img_path), "</b>"))
    })
    
    counter <- reactiveValues(
      clicknum = 0
    )
    
    img_meta <- reactiveValues(
      img = NULL,
      tempimg_path = 1, 
      first_photo = 1
    )
    
    observe({
      if((is.null(img_meta$img) || img_meta$first_photo == 1)){
        img <- jpeg::readJPEG(img_path) 
        img_meta$img <- thin2(img, thin_val)
      } 
      if(img_meta$first_photo == 1){
        z <- img_meta$img
        tempimg_path <- tempfile(fileext = ".jpg")
        jpeg::writeJPEG(z, tempimg_path)
        img_meta$tempimg_path <- tempimg_path
      }
    })
    
    
    
    
    
    click_data <- reactiveValues(
      "x" = NA,
      "y" = NA
    )
    
    observeEvent(input$image_dblclick, {
      img_meta$first_photo <- 0
      counter$clicknum <- counter$clicknum + 1
      click_data$x[counter$clicknum] <- input$image_dblclick$coords_img$y
      click_data$y[counter$clicknum] <- input$image_dblclick$coords_img$x
      
      img_meta$img <- add_point(
        img_meta$img,
        r = anchor_size_val,
        x = input$image_dblclick$coords_img$y,
        y = input$image_dblclick$coords_img$x,
        color = "blue"
      )
      
      z <- img_meta$img
      tempimg_path <- tempfile(fileext = ".jpg")
      jpeg::writeJPEG(z, tempimg_path)
      img_meta$tempimg_path <- tempimg_path
    })
    
    
    
    output$clicker <- renderPrint({
      data.frame("x" = click_data$x,
                 "y" = click_data$y)
    })
    
    
    observeEvent(input$remove_anchors, {
      img <- jpeg::readJPEG(img_path) 
      img_meta$img <- thin2(img, thin_val)
      z <- img_meta$img
      tempimg_path <- tempfile(fileext = ".jpg")
      jpeg::writeJPEG(z, tempimg_path)
      img_meta$tempimg_path <- tempimg_path
      click_data$x <- NA
      click_data$y <- NA
      img_meta$first_photo <- 1
      counter$clicknum = 0
    })
    
    
    savebutton <- reactiveValues(on = FALSE)
    observeEvent(input$save, {
      savebutton$on <- TRUE
    })
    
    observeEvent(input$reproject_button, {
      img <- jpeg::readJPEG(img_path) 
      
      pts <- list(
        list(click_data$x[1] * thin_val, click_data$y[1] * thin_val), 
        list(click_data$x[2] * thin_val, click_data$y[2] * thin_val),
        list(click_data$x[3] * thin_val, click_data$y[3] * thin_val),
        list(click_data$x[4] * thin_val, click_data$y[4] * thin_val)
      )
      
      pts <- anchor_reorder(pts)
      
      z <- reproject_grid(suppressWarnings(as.cimg(img)), 
                             init_pts = pts, 
                             dest_size = 800, 
                             qc_plot = FALSE)
      
      img_mask <- resize(as.cimg(matrix(rep(c(rep(c(1,0), 6),rep(c(0,1), 6)), 6), nrow = 12, ncol = 12)),
                         size_x = dim(z)[1], 
                         size_y = dim(z)[2])
      
      img_mask <- imappend(imlist(img_mask, img_mask, img_mask), axis = "c")
      z <- imager::imdraw(z, img_mask, opacity = 0.05)
      dim(z) <- dim(z)[-3]
      tempimg_path <- tempfile(fileext = ".jpg")
      jpeg::writeJPEG(z, tempimg_path)
      img_meta$tempimg_path <- tempimg_path
    })
    
    observeEvent(savebutton$on,{
      if(savebutton$on){
        pts <- list(
          list(click_data$x[1] * thin_val, click_data$y[1] * thin_val), 
          list(click_data$x[2] * thin_val, click_data$y[2] * thin_val),
          list(click_data$x[3] * thin_val, click_data$y[3] * thin_val),
          list(click_data$x[4] * thin_val, click_data$y[4] * thin_val)
        )
        
        out <- anchor_reorder(pts)
        attr(out, "source") <- img_path
        print(pt_list2df(out))
        assign("anchor_coord_foooooooooooooodnosdnvodnvonv", 
               value = out, 
               envir = globalenv())
        
        stopApp()
      }
    })
    
    output$plot <- renderImage({
      list(src = img_meta$tempimg_path)
    }, deleteFile = TRUE)
    
    session$onSessionEnded(function(foo) {
      stopApp()
    })
  }
  
  
  

  
  stopifnot(file.exists(img_path))
  
  global_env_variables <- c("img_path", "thin_val", "anchor_size_val",
                            "thin2","interpolate_pts","add_point")
  global_env_variables_exists <- sapply(global_env_variables,
                                        exists)
  has_conflicts <- any(global_env_variables_exists)

  if(has_conflicts){
    temp_save <- lapply(global_env_variables[global_env_variables_exists], get)
  }
  
  # Pass variables to shiny app via global environment
  anchor_size_val <<- anchor_size
  thin_val <<- thin
  img_path <<- img_path
  
  thin2 <<- function(img, thin){
    n.x <- floor(nrow(img)/thin)
    n.y <- floor(ncol(img)/thin)
    indices <- c(TRUE, rep(FALSE, (thin - 1)))
    x <- rep(indices, n.x)
    y <- rep(indices, n.y)
    img[x, y, ]
  }
  
  
  
  interpolate_pts <<- function(x,y){
    unique(
      do.call(
        "rbind",
        lapply(seq_len(length(x)-1), function(i, x, y){
          dx <- x[i+1] - x[i]
          dy <- y[i+1] - y[i]
          z <- seq(0, sqrt((dx)^2 + (dy)^2))
          theta <- atan2(dy,dx)
          
          xi <- z * cos(theta) + x[i]
          yi <- z * sin(theta) + (c(y[i]))
          return(cbind(round(xi),round(yi)))
        }, x = x, y = y)
      )
    )
  }
  
  
  add_point <<- function(img, r, x, y, color){
    color <- col2rgb(color)/255
    
    xmax <- pmin(x + r, nrow(img))
    xmin <- pmax(x - r, 1)
    ymax <- pmin(y + r, ncol(img))
    ymin <- pmax(y - r, 1)
    
    ptd <- expand.grid("xi" = (xmin:xmax), "yi" = ymin:ymax)
    v <- sqrt((ptd$xi - x)^2 + (ptd$yi - y)^2) < r
    ptd <- ptd[v,]
    v <- (ptd$xi) + nrow(img) * (ptd$yi - 1) 
    
    img[,,1][v] <- color[1,1]
    img[,,2][v] <- color[2,1]
    img[,,3][v] <- color[3,1]
    
    return(img)
  }
  
  
  
  app <- shinyApp(ui = ui_anchor_picker, 
                  server = server_anchor_picker, 
                  onStart = function(){
                    cat("\nBegin scoring...")
                    onStop(function(x){
                      cat("\nExiting scoring app...")
                    })
                  })
  shiny::runApp(app)
  
  # Remove objects created in global environment
  rm(list = global_env_variables, envir = globalenv())
  
  # Restore pre-existing global objects
  if(has_conflicts){
    for(i in seq_len(sum(global_env_variables_exists))){
      assign(eval(global_env_variables[global_env_variables_exists][i]), 
             temp_save[[i]], 
             envir = globalenv())
    }
  }
  z <- tryCatch(
    get("anchor_coord_foooooooooooooodnosdnvodnvonv", envir = globalenv()), 
    error = function(e){
      message("\nAnchor picking unsuccessful")
      return(NULL)
    }
  )
  
  tryCatch(
    rm("anchor_coord_foooooooooooooodnosdnvodnvonv", envir = globalenv()),
    warning = function(w){
      
    }
  )
  
  invisible(z)
}