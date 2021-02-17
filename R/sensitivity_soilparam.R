#' Sensitivity Analysis of soil-hydraulic parameters
#'
#' This function conducts a sensitivity analysis, comparing the effect of a max and min value of soil hydraulic input parameters
#'
#' @param project.path folder where the model input files for the current run are stored
#' @param run name of the current run, e.g. "Alfa max"
#' @param obs_nodes_all depths of the observation nodes
#' @param para list of soil-hydraulic parameters to be used in the current run
#' @param results data.frame to store the values of the soil-hydraulic parameters and the resulting RMSE and CEFF
#' @param i line in data.frame "results"
#' @return 2 figures: "Obs_points_'run'" shows timeseries of predicted soil moisture at observation nodes, "Mod_eval_run_'run'" plots observed and predicted soil moisture, with RMSE and CEFF, results
#' @export


sensitivity_soilparam <- function(project.path, obs_nodes_all, para,i, results,...){


        # write in selector.in
        write.hydraulic.para(project.path,
                             model = 0,
                             hysteresis = 0,
                             para)

        # Run hydrus model ####
        call.H1D(project.path,
                 hydrus.path = "C:/Program Files (x86)/PC-Progress/Hydrus-1D 4.xx",
                 show.output = TRUE)

        # Plot resulting soil moisture ####
        # attention: here theta_1...7 are observation nodes i.e. 5,10,20,30,40,60,100 cm!!!
        pred <- read.obs_node(project.path,
                              out.file = "Obs_Node.out",
                              obs.output = "theta",
                              obs.nodes = obs_nodes_all)
        png(paste("Figures/Obs_points_", run,".png"), height=5, width=7, units = "in", res=500,)
        plot(pred$Time, pred$theta_5, type = "l", col = 1, ylim = c(0,1),
             xlab = "Day", ylab = "Soil moisture (m^3/m^3)")
        lines(pred$Time, pred$theta_10, type = "l", col = 2)
        lines(pred$Time, pred$theta_20, type = "l", col = 3)
        lines(pred$Time, pred$theta_30, type = "l", col = 4)
        lines(pred$Time, pred$theta_40, type = "l", col = 5)
        lines(pred$Time, pred$theta_60, type = "l", col = 6)
        lines(pred$Time, pred$theta_100, type = "l", col = 7)

        legend("bottomright", legend = c("5cm", "10cm", "20cm", "30cm", "40cm", "60cm", "100cm"),
               col = c(1:7), lty = 1)
        title(run)
        dev.off()

        # Calculate quality measures ####
        # RMSE and Nash-Sutcliffe CEFF simulation/observation (obs nodes 2,3,4,5)
        # attention: here theta_1...4 are points observed by PR2-probe i.e. 10,20,30,40 cm!!!
        pr <- read.table("C:/Dateien/Arbeit/WiMi_2020/Bewamo/process_rhinluch/rhi_pr2_theta.txt", sep = "\t")
        pr <- aggregate(x = pr[,7:10], by = list(as.Date(pr$datetime)), FUN = mean)
        colnames(pr)[1] <- "date"
        pr <- pr[1:(nrow(pr)-1),]
        pr$date <- 1:nrow(pr)
        data <- gather(pr, key = "probe", value = "obs","theta_1", "theta_2", "theta_3", "theta_4")

        predtidy <- pred[pred$Time %in% 1:nrow(pr),c(1,3:6)]
        colnames(predtidy) <- c("date", "theta_1","theta_2", "theta_3", "theta_4")
        predtidy <- gather(predtidy, key = "probe", value =  "pred","theta_1", "theta_2", "theta_3", "theta_4")

        data$pred <- predtidy$pred

        png(paste("Figures/Mod_eval_", run,".png"), height=5, width=7, units = "in", res=500)
        mod_eval <- model.eval(data,
                               var.name = "theta",
                               print.stats = T,
                               plot = T)
        legend("topleft", legend = c(names(mod_eval), round(mod_eval, digits = 2)), ncol = 2, bty = "n")
        title(run)
        dev.off()

        # mod_eval <- model.eval(data, var.name = "theta", print.stats = T, plot = F)

        results$thr[i] <- para$thr
        results$ths[i] <- para$ths
        results$Alfa[i] <- para$Alfa
        results$n[i] <- para$n
        results$Ks[i] <- para$Ks
        results$l[i] <- para$l
        results$RMSE[i] <- mod_eval["RMSE"]
        results$CEFF[i] <- mod_eval["CEFF"]

        return(results)

}
