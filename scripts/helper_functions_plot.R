#'############################################################################'#
#'################## HELPER FUNCTIONS FOR DATA VISUALIZATION #################'#
#'############################################################################'#

library(plotly)
library(shiny)
library(mapview)
library(sp)
library(sf)
# library(htmlwidgets)
library(maptools)
library(ggforce)
library(ggplot2)
library(patchwork)

### GRAPHICAL PARAMETERS for MAPS IN mapview -----------------------------------
map.type = "Esri.WorldTopoMap"
col_mesh = rgb(60, 60, 60, max = 255)
col_bdy = rgb(128, 128, 128, max = 255)
col_nodes = rgb(0, 0, 0, max = 255)
col_locations = rgb(190, 0, 0, max = 255)

### PLOT MESH 2D USING plotly --------------------------------------------------
plot.mesh.2D <- function(mesh, ...){
  
  # Plot 
  p = plot_ly() %>% layout(scene = list(
    aspectratio = list(
      x = 1,
      y = 1
    )),
    xaxis = list(
      scaleanchor = "y",
      title = "",
      showgrid = F,
      zeroline = F,
      showticklabels = F,
      ticks = ""),
    yaxis = list(
      scaleanchor = "x",
      title = "",
      showgrid = F,
      zeroline = F,
      showticklabels = F,
      ticks = ""),
    margin = list(
      b = 0,
      l = 0,
      r = 0,
      t = 0
    )) %>%
    add_segments(x = mesh$nodes[mesh$edges[,1],1],
                 y = mesh$nodes[mesh$edges[,1],2],
                 xend = mesh$nodes[mesh$edges[,2],1],
                 yend = mesh$nodes[mesh$edges[,2],2],
                 color = I("gray"),
                 showlegend = F) %>%
    add_markers(x = mesh$nodes[,1],
                y = mesh$nodes[,2],
                color = I("gray2"),
                marker = list(size = 5),
                hoverinfo = "text",
                text = paste("</br><b> Coordinates:", round(mesh$nodes[,1],2),
                             round(mesh$nodes[,2],2)),
                showlegend = F,
                visible = T)
  
  return(p)
  
}


### PLOT SAMPLE 2D USING plotly ------------------------------------------------
plot.sample.2D <- function(data, mesh, ...){
  
  # Plot
  p = plot_ly(type = "scatter", mode = "markers") %>%
    layout(scene = list(
      aspectratio = list(x = 1, y = 1)
    ),
    xaxis = list(
      scaleanchor = "y",
      title = "",
      showgrid = F,
      zeroline = F,
      showticklabels = F,
      ticks = ""),
    yaxis = list(
      scaleanchor = "x",
      title = "",
      showgrid = F,
      zeroline = F,
      showticklabels = F,
      ticks = ""),
    margin = list(
      b = 0,
      l = 0,
      r = 0,
      t = 0
    )) %>%
    add_segments(x = mesh$nodes[mesh$edges[,1],1],
                 y = mesh$nodes[mesh$edges[,1],2],
                 xend = mesh$nodes[mesh$edges[,2],1],
                 yend = mesh$nodes[mesh$edges[,2],2],
                 color = I("gray"),
                 showlegend = F) %>%
    add_markers(x = mesh$nodes[,1],
                y = mesh$nodes[,2],
                color = I("gray2"),
                marker = list(size = 2.5),
                opacity = 1,
                hoverinfo = "text",
                text = paste("</br><b> Coordinates:", round(mesh$nodes[,1],2),
                             round(mesh$nodes[,2],2)),
                showlegend = F,
                visible = T) %>%
    add_markers(x = data[,1],
                y = data[,2],
                color = I("red3"),
                marker = list(size = 5),
                opacity = 1,
                hoverinfo = "text",
                text = paste("</br><b> Coordinates:", round(data[,1],2),
                             round(data[,2],2)),
                showlegend = F,
                visible = T)
  
  return(p)
  
}


## PLOT DENSITY 2D USING image -------------------------------------------------
plot.density.2D <- function(X, Y, Z, min_range = NULL, max_range = NULL, boundary = NULL, colorscale = NULL, ...){
  
  if(is.null(min_range)) {min_range = min(Z, na.rm = TRUE)}
  if(is.null(max_range)) {max_range = max(Z, na.rm = TRUE)}
  if(is.null(colorscale)) {colorscale = "heat.colors"}
  color_palette <- match.fun(colorscale)
  color_palette <- color_palette(100)
  
  image2D(x = X, y = Y, z = matrix(Z, nrow = length(X), ncol = length(Y), byrow = FALSE),
          col = color_palette, zlim = c(min_range, max_range),
          asp = 1, xlab = '', ylab = '', contour = list(nlevels = 10, drawlabels = FALSE),
          axes = FALSE, frame.plot = FALSE )
  
  if(!is.null(boundary)){
    points(boundary, type = 'l', lwd = 3)
  }
  
  p = recordPlot()
  
  return(p)
  
}


### CONVERT A mesh INTO A sfc OBJECT -------------------------------------------
st_as_sfc.mesh.2D <- function(mesh, crs = NULL, ...){
  
  polygon_list <- apply(mesh$triangles, MARGIN = 1, FUN = function(elem){
    st_cast(st_linestring(mesh$nodes[elem,]), to ="POLYGON")
  })
  
  if(is.null(crs)){crs <- NA_crs_}
  mesh_sf <- st_sfc(polygon_list, crs = crs)
  
  return(mesh_sf)
  
}


### PLOT MESH 2D ON MAP USING mapview ------------------------------------------
plot.mesh.2D.map <- function(mesh, ...){
  
  # From British National Grid CRS (EPSG:27700) to latitude/longitude CRS (EPSG:4326 - WSG 84)
  mesh$nodes = mesh$nodes * 1000
  mesh_df = as.data.frame(mesh$nodes[,2:1])
  colnames(mesh_df) = c("Easting", "Northing")
  mesh_sp = mesh_df
  coordinates(mesh_sp) = ~Northing+Easting
  proj4string(mesh_sp) = CRS("+init=epsg:27700")
  mesh_sp = spTransform(mesh_sp, CRS("+init=epsg:4326"))
  
  mesh$nodes = mesh_sp@coords
  mesh_df = data.frame(lon = mesh$nodes[,1], lat = mesh$nodes[,2])
  
  mesh_sf = st_as_sfc.mesh.2D(mesh, crs = 4326)
  mesh_sf = st_as_sf(mesh_sf)
  
  map = mapview(mesh_sf, legend = FALSE, map.type = map.type,
                color = col_mesh, col.regions = col_mesh,
                layer.name = "mesh", lwd = 1, alpha.regions = 0.25) # +
  # mapview(st_as_sf(mesh_df, coords = c("lon", "lat"), crs = 4326),
  #         legend = FALSE, col.region = col_nodes, layer.name = "mesh-nodes",
  #         alpha.regions = 1, cex = 1.5)
  
  return(map)
  
}


### PLOT SAMPLE 2D ON MAP USING mapview ----------------------------------------
plot.sample.2D.map <- function(data, mesh, ...){
  
  # From British National Grid CRS (EPSG:27700) to latitude/longitude CRS (EPSG:4326 - WSG 84)
  mesh$nodes = mesh$nodes * 1000
  mesh_df = as.data.frame(mesh$nodes[,2:1])
  colnames(mesh_df) = c("Easting", "Northing")
  mesh_sp = mesh_df
  coordinates(mesh_sp) = ~Northing+Easting
  proj4string(mesh_sp) = CRS("+init=epsg:27700")
  mesh_sp = spTransform(mesh_sp, CRS("+init=epsg:4326"))
  mesh$nodes = mesh_sp@coords
  mesh_df = data.frame(lon = mesh$nodes[,1], lat = mesh$nodes[,2])
  mesh_sf = st_as_sfc.mesh.2D(mesh, crs = 4326)
  mesh_sf = st_as_sf(mesh_sf)
  
  data_df = data.frame(Easting = data[,2] * 1000, Northing = data[,1] * 1000)
  data_sp = data_df
  coordinates(data_sp) = ~Northing+Easting
  proj4string(data_sp) = CRS("+init=epsg:27700")
  data_sp = spTransform(data_sp, CRS("+init=epsg:4326"))
  data_df = data_sp@coords
  data_df = as.data.frame(data_df)
  names(data_df) = c("lon", "lat")
  
  map = mapview(mesh_sf, legend = FALSE, map.type = map.type,
                color = col_mesh, col.regions = col_mesh,
                layer.name = "mesh", lwd = 1, alpha.regions = 0.25) +
    # mapview(st_as_sf(mesh_df, coords = c("lon", "lat"), crs = 4326),
    #         legend = FALSE, col.region = col_nodes, layer.name = "mesh-nodes",
    #         alpha.regions = 1, cex = 1.5) +
    mapview(st_as_sf(data_df, coords = c("lon", "lat"), crs = 4326),
            legend = FALSE, col.region = col_locations, layer.name = "data-locations",
            alpha.regions = 1, cex = 2)
  
  return(map)
  
}


### PLOT DENSITY 2D USING mapview ----------------------------------------------
plot.density.2D.map <- function(coeff, mesh, colorscale = NULL, ...){
  
  coeff <- apply(mesh$triangles, MARGIN = 1, FUN = function(edge){
    mean(coeff[edge])
  })
  
  if(is.null(colorscale)) {colorscale = "heat.colors"}
  color_palette <- match.fun(colorscale)
  color_palette <- color_palette(100)
  
  # From British National Grid CRS (EPSG:27700) to latitude/longitude CRS (EPSG:4326 - WSG 84)
  mesh$nodes = mesh$nodes * 1000
  mesh_df = as.data.frame(mesh$nodes[,2:1])
  colnames(mesh_df) = c("Easting", "Northing")
  mesh_sp = mesh_df
  coordinates(mesh_sp) = ~Northing+Easting
  proj4string(mesh_sp) = CRS("+init=epsg:27700")
  mesh_sp = spTransform(mesh_sp, CRS("+init=epsg:4326"))
  mesh$nodes = mesh_sp@coords
  mesh_df = data.frame(lon = mesh$nodes[,1], lat = mesh$nodes[,2])
  mesh_sf = st_as_sfc.mesh.2D(mesh, crs = 4326)
  
  U = st_sf(data.frame(coeff = coeff), geometry = mesh_sf)
  
  map = mapview(U, legend = TRUE, map.type = map.type,
                color = col_mesh, col.regions = color_palette,
                layer.name = "density", alpha.regions = 1)
  
  return(map)
  
}


### PLOT MESH 2.5D USING rgl ---------------------------------------------------
plot.mesh.2.5D <- function(mesh, ...){
  
  coeff = rep(0, nrow(mesh$nodes))
  
  plot.density.2.5D(coeff, mesh, ...)
  
  rglwidget()
}


### PLOT SAMPLE 2.5D USING rgl -------------------------------------------------
plot.sample.2.5D <- function(data, mesh, radius = NULL, ...){
  
  if(is.null(radius)){radius = 0.015}
  
  plot.mesh.2.5D(mesh = mesh, ...)
  
  rgl.spheres(x = data[,1], y = data[,2], z = data[,3], radius = radius,
              color = "red3")
  
  rglwidget()
}


## FEM PLOT (Colormap with Range [m,M]) ----------------------------------------
plot.density.2.5D <- function(coeff, mesh, M = NULL, m = NULL, colorscale = NULL, world = FALSE, ...){
  
  FEM = FEM(coeff = coeff, create.FEM.basis(mesh = mesh))
  
  if (is.null(m)) {m = min(FEM$coeff)}
  if (is.null(M)) {M = max(FEM$coeff)}
  triangles = c(t(FEM$FEMbasis$mesh$triangles))
  ntriangles = nrow(FEM$FEMbasis$mesh$triangles)
  order = FEM$FEMbasis$mesh$order
  nodes = FEM$FEMbasis$mesh$nodes
  edges = matrix(rep(0, 6*ntriangles), ncol = 2)
  for(i in 0:(ntriangles-1)){
    edges[3*i+1,] = c(triangles[3*order*i+1], triangles[3*order*i+2])
    edges[3*i+2,] = c(triangles[3*order*i+1], triangles[3*order*i+3])
    edges[3*i+3,] = c(triangles[3*order*i+2], triangles[3*order*i+3])
  }
  edges = edges[!duplicated(edges),]
  edges <- as.vector(t(edges))
  
  coeff = FEM$coeff
  FEMbasis = FEM$FEMbasis
  mesh = FEMbasis$mesh
  
  if(is.null(colorscale)){colorscale = "heat.colors"}
  color_palette <- match.fun(colorscale)
  color_palette <- color_palette(100)
  ncolor = length(color_palette)
  palette(color_palette)
  
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d()
    rgl.pop("lights") 
    light3d(specular = "black")
    
    diffrange = M - m
    
    col = coeff[triangles,isurf]
    col = (col - min(coeff[,isurf]))/diffrange*(ncolor-1)+1
    
    rgl.triangles(x = nodes[triangles,1], y = nodes[triangles,2],
                  z = nodes[triangles,3],
                  color = col,...)
    rgl.lines(x = nodes[edges,1], y = nodes[edges,2],
              z = nodes[edges,3],
              color = "gray40",...)
    aspect3d("iso")
    
    if (nsurf > 1 && isurf<nsurf)
    {readline("Press a button for the next plot...")}
  }
  
  # Import world coastlines
  xyz <- shapefile_lines()
  
  # Add world coastlines to the plot
  if(world){rgl.linestrips(x = xyz[,1], y = xyz[,2], z = xyz[,3], col = "black", lwd = 3)}
  
  rglwidget()
}


## WORLD SHAPEFILE LINES -------------------------------------------------------
shapefile_lines <- function() {
  
  shp <- readShapeSpatial(fn = "data/coastlines/ne_110m_coastline.shp")
  
  # Combine lines in a single matrix
  mat <- do.call(rbind, sapply(1:length(shp@lines), function(i) rbind(shp@lines[[i]]@Lines[[1]]@coords,c(NA,NA))))
  
  # Convert spherical to cartesian
  xyz <- rgl.sph2car(mat[,2], mat[,1], radius = 1)
  
  return (xyz)
}

# From spherical to Cartesian coordinates
rgl.sph2car <- function(lat = 0, lon = 0, radius = 1, deg = T, precise = T) {
  if(deg) {
    if(precise){
      lat <- lat/180
      lon <- lon/180
      x <- radius*cospi(lat)*cospi(lon)
      y <- radius*cospi(lat)*sinpi(lon)
      z <- radius*sinpi(lat)
    }else{
      lat <- lat*pi/180
      lon <- lon*pi/180
      x <- radius*cos(lat)*cos(lon)
      y <- radius*cos(lat)*sin(lon)
      z <- radius*sin(lat)
    }
  }
  return(matrix(c(x,y,z), nrow = length(x), ncol = 3, dimnames = list(NULL, c("x","y","z"))))
}


### PLOT MESH 1.5D USING ggplot ------------------------------------------------
plot.mesh.1.5D <- function(mesh, ...){
  
  num_edges = dim(mesh$edges)[1]
  
  x = vector(mode = "double", length = 2*num_edges)
  y = vector(mode = "double", length = 2*num_edges)
  grp.nodes = vector(mode = "integer", length = 2*num_edges)
  
  for(e in 1:num_edges){
    x[(2*(e-1)+1):(2*(e-1)+2)] = c(mesh$nodes[mesh$edges[e,1], 1], mesh$nodes[mesh$edges[e,2],1])
    y[(2*(e-1)+1):(2*(e-1)+2)] = c(mesh$nodes[mesh$edges[e,1], 2], mesh$nodes[mesh$edges[e,2],2])
    grp.nodes[(2*(e-1)+1):(2*(e-1)+2)] = rep(e, times = 2)
  }
  
  data_plot = data.frame(x = x, y = y, group = grp.nodes)
  
  margin_size = 0.1
  
  ggplot(data_plot) +
    geom_link2(aes(x = x, y = y, group = grp.nodes),
               lineend = 'round', n = 1, col = "gray50", ...) +
    labs(x = "", y = "", color = "", title = "") +  
    coord_fixed(ratio = 1) + theme_void() +
    theme(plot.margin = unit(c(margin_size, margin_size, margin_size, margin_size), "cm")) +
    geom_point(data = data.frame(x = mesh$nodes[,1], y = mesh$nodes[,2]),
               aes(x = x, y = y), color = "black", size = 2)
  
}


### PLOT SAMPLE 1.5D USING ggplot ----------------------------------------------
plot.sample.1.5D <- function(data, mesh, ...){
  
  num_edges = dim(mesh$edges)[1]
  
  x = vector(mode = "double", length = 2*num_edges)
  y = vector(mode = "double", length = 2*num_edges)
  grp.nodes = vector(mode = "integer", length = 2*num_edges)
  
  for(e in 1:num_edges){
    x[(2*(e-1)+1):(2*(e-1)+2)] = c(mesh$nodes[mesh$edges[e,1], 1], mesh$nodes[mesh$edges[e,2],1])
    y[(2*(e-1)+1):(2*(e-1)+2)] = c(mesh$nodes[mesh$edges[e,1], 2], mesh$nodes[mesh$edges[e,2],2])
    grp.nodes[(2*(e-1)+1):(2*(e-1)+2)] = rep(e, times = 2)
  }
  
  data_plot = data.frame(x = x, y = y, group = grp.nodes)
  
  margin_size = 0.1
  
  ggplot(data_plot) +
    geom_link2(aes(x = x, y = y, group = grp.nodes),
               lineend = 'round', n = 1, col = "gray50", ...) +
    labs(x = "", y = "", color = "", title = "") +  
    coord_fixed(ratio = 1) + theme_void() +
    theme(plot.margin = unit(c(margin_size, margin_size, margin_size, margin_size), "cm")) +
    geom_point(data = data.frame(x = data[,1], y = data[,2]),
               aes(x = x, y = y), color = "red3", size = 2)
  
}


### PLOT DENSITY 1.5D USING ggplot ---------------------------------------------
plot.density.1.5D <- function(coeff, mesh, m = NULL, M = NULL,
                              colorscale = NULL, ...){
  
  FEM = FEM(coeff = coeff, create.FEM.basis(mesh = mesh))
  num_edges = dim(mesh$edges)[1]
  
  if (is.null(m)) {m = min(FEM$coeff, na.rm = TRUE)}
  if (is.null(M)) {M = max(FEM$coeff, na.rm = TRUE)}
  if(is.null(colorscale)){colorscale = "heat.colors"}
  color_palette <- match.fun(colorscale)
  color_palette <- color_palette(100)
  
  x = vector(mode = "double", length = 2*num_edges)
  y = vector(mode = "double", length = 2*num_edges)
  coeff = vector(mode = "double", length = 2*num_edges)
  grp.nodes = vector(mode = "integer", length = 2*num_edges)
  
  for(e in 1:num_edges){
    x[(2*(e-1)+1):(2*(e-1)+2)] = c(mesh$nodes[mesh$edges[e,1], 1], mesh$nodes[mesh$edges[e,2],1])
    y[(2*(e-1)+1):(2*(e-1)+2)] =  c(mesh$nodes[mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2],2])
    coeff[(2*(e-1)+1):(2*(e-1)+2)] = c(FEM$coeff[mesh$edges[e,1]], FEM$coeff[mesh$edges[e,2]])  
    grp.nodes[(2*(e-1)+1):(2*(e-1)+2)] = rep(e,times=2)
  }
  
  data_plot = data.frame(x = x, y = y, coeff = coeff, grp.nodes)
  
  margin_size = 0.1
  plot_height = 2
  
  plot = ggplot(data_plot) +
    geom_link2(aes(x = x, y = y, colour = coeff, group = grp.nodes),
               lineend = 'round', n = 10, ...) +
    labs(x = "", y = "", color = "", title = "") +  
    coord_fixed(ratio = 1) + theme_void() +
    theme(plot.margin = unit(c(margin_size, margin_size, margin_size, margin_size), "cm")) +
    scale_color_gradientn(colors = color_palette, limits = c(m, M))
  
  plot = plot + theme(legend.position = "right",
                      legend.direction = "vertical",
                      legend.box = "vertical",
                      legend.key.width = unit(0.25, "in"),
                      legend.key.height = unit(1, "in"))
  
  return(plot)
  
}


### PLOT MESH 1.5D USING mapview -----------------------------------------------
plot.mesh.1.5D.map <- function(mesh, bdy = NULL, ...){
  
  mesh_df = data.frame(lon = mesh$nodes[,1], lat = mesh$nodes[,2])
  mesh_df = st_as_sf(mesh_df, coords = c("lon", "lat"), crs = 4326)
  mesh_linnet = as.linnet(mesh)
  mesh_sfnetwork = as_sfnetwork(mesh_linnet, directed = FALSE, edges_as_lines = TRUE)
  st_crs(mesh_sfnetwork) = 4326
  mesh_sf = st_as_sf(mesh_sfnetwork, "edges")
  
  plot = mapview(mesh_sf, legend = FALSE, map.type = map.type, color = col_mesh,
                 layer.name = "road-network", lwd = 2, alpha.regions = 1) +
    mapview(mesh_df, legend = FALSE, col.region = col_nodes, layer.name = "mesh-nodes",
            alpha.regions = 1, cex = 1)
  
  if(!is.null(bdy)){
    plot = plot + mapview(st_sf(data.frame(x = 1), geoemtry = bdy), legend = FALSE,
                          col.region = col_bdy, alpha.regions = 0.3, layer.name = "boundary", lwd = 0.3) 
  }
  
  return(plot)
  
}


### PLOT SAMPLE 1.5D USING mapview ---------------------------------------------
plot.sample.1.5D.map <- function(mesh, data, bdy = NULL, ...){
  
  data = st_as_sf(data, coords = c("lon", "lat"), crs = 4326)
  mesh_linnet = as.linnet(mesh)
  mesh_sfnetwork = as_sfnetwork(mesh_linnet, directed = FALSE, edges_as_lines = TRUE)
  st_crs(mesh_sfnetwork) = 4326
  mesh_sf = st_as_sf(mesh_sfnetwork, "edges")
  
  plot = mapview(mesh_sf, legend = FALSE, map.type = map.type, color = col_mesh,
                 layer.name = "road-network", lwd = 2, alpha.regions = 1) +
    mapview(data, legend = FALSE, col.region = col_locations, layer.name = "data",
            alpha.regions = 1, cex = 2.5)
  
  if(!is.null(bdy)){
    plot = plot + mapview(st_sf(data.frame(x = 1), geoemtry = bdy), legend = FALSE,
                          col.region = col_bdy, alpha.regions = 0.3, layer.name = "boundary", lwd = 0.3) 
  }
  
  return(plot)
  
}


### PLOT DENSITY 1.5D USING mapview --------------------------------------------
plot.density.1.5D.map <- function(coeff, mesh, bdy = NULL, colorscale = NULL, ...){
  
  mesh_linnet = as.linnet(mesh)
  mesh_sfnetwork = as_sfnetwork(mesh_linnet, directed = FALSE, edges_as_lines = TRUE)
  st_crs(mesh_sfnetwork) = 4326
  mesh_sf = st_as_sf(mesh_sfnetwork, "edges")
  
  coeff_edges = vector(mode = "numeric", length = nrow(mesh$edges))
  for(e in 1:nrow(mesh$edges)){
    coeff_edges[e] = mean(c(coeff[mesh$edges[e,1]], coeff[mesh$edges[e,2]]), na.rm = TRUE)
  }
  
  U = st_sf(data.frame(coeff = coeff_edges), geometry = st_as_sfc(mesh_sf))
  
  if(is.null(colorscale)) {colorscale = "heat.colors"}
  color_palette <- match.fun(colorscale)
  color_palette <- color_palette(100)
  
  map = mapview(U, legend = TRUE, map.type = map.type,
                color = color_palette, col.regions = color_palette,
                layer.name = "density", alpha.regions = 1, lwd = 2)
  
  if(!is.null(bdy)){
    map = map + mapview(st_sf(data.frame(x = 1), geoemtry = bdy), legend = FALSE,
                        col.region = col_bdy, alpha.regions = 0.3, layer.name = "boundary", lwd = 0.3) 
  }
  
  return(map)
  
}