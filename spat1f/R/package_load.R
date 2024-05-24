
load_all2 <- function (path = ".", reset = TRUE, recompile = FALSE, export_all = TRUE, 
         helpers = TRUE, quiet = FALSE, ...) {
  if (inherits(path, "package")) {
    path <- path$path
  }
  devtools:::save_all()
  
  .load_all_inernal <- function (path = ".", reset = TRUE, 
                                 compile = NA, attach = TRUE, 
                                 export_all = TRUE, export_imports = export_all, 
                                 helpers = export_all, 
                                 attach_testthat = pkgload:::uses_testthat(path), 
                                 quiet = NULL, recompile = FALSE, 
                                 warn_conflicts = TRUE) {
    path <- pkgload:::pkg_path(path)
    package <- pkgload:::pkg_name(path)
    description <- pkgload:::pkg_desc(path)
    withr::local_envvar(c(DEVTOOLS_LOAD = package))
    quiet <- pkgload:::load_all_quiet(quiet, "load_all")
    if (!quiet) {
      cli::cli_inform(c(i = "Loading {.pkg {package}}"))
    }
    if (package == "compiler") {
      oldEnabled <- compiler::enableJIT(0)
      on.exit(compiler::enableJIT(oldEnabled), TRUE)
    }
    if (missing(compile) && !missing(recompile)) {
      compile <- if (isTRUE(recompile)) 
        TRUE
      else NA
    }
    if (isTRUE(compile)) {
      pkgbuild::clean_dll(path)
      pkgbuild::compile_dll(path, quiet = quiet)
    }
    else if (identical(compile, NA)) {
      pkgbuild::compile_dll(path, quiet = quiet)
    }
    else if (identical(compile, FALSE)) {
    }
    else {
      cli::cli_abort("{.arg compile} must be a logical vector of length 1.")
    }
    old_methods <- list()
    if (reset) {
      pkgload:::clear_cache()
      if (pkgload:::is_loaded(package)) {
        pkgload:::patch_colon(package)
        methods_env <- pkgload:::ns_s3_methods(package)
        pkgload:::unregister(package)
        old_methods <- as.list(methods_env)
        old_methods <- Filter(function(x) pkgload:::is_foreign_method(x, 
                                                                      package), old_methods)
      }
    }
    if (pkgload:::is_loaded(package)) {
      rlang::env_unlock(pkgload:::ns_env(package))
    }
    else {
      pkgload:::create_ns_env(path)
    }
    out <- list(env = pkgload:::ns_env(package))
    pkgload:::load_depends(path, quiet = quiet)
    pkgload:::load_imports(path)
    pkgload:::insert_imports_shims(package)
    out$data <- pkgload:::load_data(path)
    out$code <- pkgload:::load_code(path, quiet = quiet)
    pkgload:::register_s3(path)
    if (identical(compile, FALSE)) {
      out$dll <- pkgload:::try_load_dll(path)
    }
    else {
      out$dll <- .load_dll2(path)
    }
    if (isTRUE(attach_testthat) && package != "testthat") {
      ("base" %:::% "library")("testthat", warn.conflicts = FALSE)
    }
    pkgload:::load_po(package, path)
    pkgload:::run_pkg_hook(package, "load")
    pkgload:::setup_ns_exports(path)
    pkgload:::run_ns_load_actions(package)
    ns <- pkgload:::ns_env(package)
    lockEnvironment(ns)
    for (nm in names(ns)) {
      lockBinding(nm, ns)
    }
    pkgload:::run_user_hook(package, "load")
    if (attach) {
      pkgload:::setup_pkg_env(package)
    }
    rlang::env_bind(pkgload:::ns_s3_methods(package), !!!old_methods)
    if (attach) {
      pkgload:::run_pkg_hook(package, "attach")
      pkgload:::run_user_hook(package, "attach")
      pkgload:::populate_pkg_env(package, path, export_all, export_imports, 
                                 helpers)
    }
    pkgload:::insert_global_shims()
    if (isTRUE(warn_conflicts)) {
      pkgload:::warn_if_conflicts(package, out$env, globalenv())
    }
    invisible(out)
  }
  
  .load_dll2 <- function(path = "."){
    package <- pkgload:::pkg_name(path)
    env <- pkgload:::ns_env(package)
    nsInfo <- pkgload:::parse_ns_file(path)
    dlls <- list()
    dynLibs <- nsInfo$dynlibs
    nativeRoutines <- list()
    for (i in seq_along(dynLibs)) {
      lib <- dynLibs[i]
      dlls[[lib]] <- .library.dynam3(path, lib)
      routines <- pkgload:::assignNativeRoutines(dlls[[lib]], lib, env, 
                                                 nsInfo$nativeRoutines[[lib]])
      nativeRoutines[[lib]] <- routines
      if (!is.null(names(nsInfo$dynlibs)) && nzchar(names(nsInfo$dynlibs)[i])) 
        env[[names(nsInfo$dynlibs)[i]]] <- dlls[[lib]]
      setNamespaceInfo(env, "DLLs", dlls)
    }
    pkgload:::addNamespaceDynLibs(env, nsInfo$dynlibs)
    dll_path <- dlls[[package]][["path"]]
    if (!is_null(dll_path)) {
      rlang::new_weakref(env, finalizer = pkgload:::ns_finalizer(dll_path))
    }
    invisible(dlls)
  }
  
  .library.dynam3 <- function(path = ".", lib = ""){
    path <- pkgload:::pkg_path(path)
    dyn_ext <- .Platform$dynlib.ext
    dllname <- paste(lib, dyn_ext, sep = "")
    dllfile <- pkgload:::package_file("src", dllname, path = path)
    if (!file.exists(dllfile)) {
      return(invisible())
    }
    dllinfo <- dyn.load(dllfile)
    .dynLibs(c(.dynLibs(), list(dllinfo)))
    return(dllinfo)
  }
  
  
  .load_all_inernal(path = path, reset = reset, recompile = recompile, 
                    export_all = export_all, helpers = helpers, quiet = quiet, 
                    ...)
}


