import ClimaAtmos: time_from_filename
import ClimaCore: Geometry, Spaces, Fields, InputOutput
import ClimaComms
import CairoMakie: Makie
import Statistics: mean


"""
Generic Functions which allow other plotting functions to work.
Generally file_path or folder_path, and the time interval of the saved simulation 
files are necessary to use the later plotting function
"""
function read_hdf5_file(file_path)
    reader =
        InputOutput.HDF5Reader(file_path, ClimaComms.SingletonCommsContext())
    diagnostics = InputOutput.read_field(reader, "diagnostics")
    close(reader)
    return time_from_filename(file_path), diagnostics
end

function create_time_vals(interval, max_time)
    num_intervals = Int(max_time / interval) + 1
    time_vals = Array{Float64}(undef, num_intervals)
    time_vals[1] = 0.0
    for i in 1:(num_intervals - 1)
        time_vals[i + 1] = i * Float64(interval)
    end
    return time_vals
end

function create_day_file_paths(interval, max_day)
    time_vals = create_time_vals(interval, max_day)
    day_vals = Array{String}(undef, length(time_vals))
    count = 1
    for time in time_vals
        day_vals[count] = "/day" * string(time_vals[count]) * ".hdf5"
        count += 1
    end
    return day_vals
end


"""
Examples of some plotting functions that are possible with output data,
Note: all functions plot profiles with respect to either altitude or pressure
as that is how we extract and visualize relevant data
"""

# Contour function plot of temperature over time
function plot_temp(folder_paths, graph_titles, max_temp, figureName)
    figure = Makie.Figure()
    Makie.empty!(figure)
    days_file_paths = days_file_paths = create_day_file_paths(40, 200)
    time_vals = Float32[0.0, 40.0, 80.0, 120.0, 160.0, 200.0]
    _, diagnostics_example =
        read_hdf5_file(folder_paths[1] * days_file_paths[1])
    z = vec(
        parent(
            Fields.coordinate_field(axes(diagnostics_example.temperature)).z,
        ),
    )
    plots_per_row = Int(length(folder_paths) / 2)
    graph_num = 1
    for folder_path in folder_paths
        time_count = 1
        temp = Matrix{Float32}(undef, length(days_file_paths), length(z))
        for day in days_file_paths
            file_path = folder_path * day
            _, diagnostics = read_hdf5_file(file_path)
            temp[time_count, :] = vec(parent(diagnostics.temperature))
            time_count += 1
        end
        axis = Makie.Axis(
            figure[
                trunc(Int, (graph_num - 1) / plots_per_row),
                (graph_num - 1) % plots_per_row,
            ],
            xlabel = "Time [days]",
            ylabel = "Elevation [m]",
            title = graph_titles[graph_num],
            height = 160,
            width = 165 * 3 / plots_per_row,
            xlabelsize = 14.0,
            ylabelsize = 14.0,
            xminorticksize = 2.0,
            yminorticksize = 2.0,
            xticklabelsize = 10.0,
            yticklabelsize = 10.0,
        )
        Makie.contourf!(
            axis,
            time_vals,
            z,
            temp,
            levels = 30,
            colormap = :viridis,
            extendlow = :red,
        )
        graph_num += 1
    end
    Makie.Colorbar(
        figure[2, :],
        limits = (200, max_temp),
        colormap = :viridis,
        label = "Temperature [K]",
        vertical = false,
    )
    Makie.display(figure)
    Makie.save(figureName, figure)
end

# plots the final temperature profile at the end of a simulation
function plot_final_temp_state(file_path)
    figure = Makie.Figure()
    Makie.empty!(figure)
    time, diagnostics = read_hdf5_file(file_path)
    z = vec(parent(Fields.coordinate_field(axes(diagnostics.temperature)).z))
    temp = vec(parent(diagnostics.temperature))
    axis = Makie.Axis(
        figure[1, 1],
        title = "Final Temperature Profile",
        xlabel = "Temperature [K]",
        ylabel = "Height [m]",
    )
    Makie.lines!(axis, temp, z) 
    Makie.scatter!(axis, temp[1] + 3, z[1])
    Makie.display(figure)
    Makie.save("example_temp_profile.png", figure)
end

# plots shortwave flux contour plot over time
function plot_sw_flux_down_time(
    folder_paths,
    graph_titles,
    figureName,
    low_lim,
    up_lim,
)
    figure = Makie.Figure()
    Makie.empty!(figure)
    days_file_paths = create_day_file_paths(40, 200)
    time_vals = Float32[0.0, 40.0, 80.0, 120.0, 160.0, 200.0]
    _, diagnostics_example =
        read_hdf5_file(folder_paths[1] * days_file_paths[1])
    _, z... = vec(
        parent(
            Fields.coordinate_field(axes(diagnostics_example.sw_flux_down)).z,
        ),
    )
    plots_per_row = Int(length(folder_paths) / 2)
    graph_num = 1
    for folder_path in folder_paths
        time_count = 1
        swd = Matrix{Float32}(undef, length(days_file_paths), length(z))
        for day in days_file_paths
            file_path = folder_path * day
            time, diagnostics = read_hdf5_file(file_path)
            _, swd[time_count, :]... = vec(parent(diagnostics.sw_flux_down))
            time_count += 1
        end
        axis = Makie.Axis(
            figure[
                trunc(Int, (graph_num - 1) / plots_per_row),
                (graph_num - 1) % plots_per_row,
            ],
            xlabel = "Time [days]",
            ylabel = "Elevation [m]",
            title = graph_titles[graph_num],
            height = 160,
            width = 165 * 3 / plots_per_row,
            xlabelsize = 14.0,
            ylabelsize = 14.0,
            xminorticksize = 2.0,
            yminorticksize = 2.0,
            xticklabelsize = 10.0,
            yticklabelsize = 10.0,
        )
        Makie.contourf!(
            axis,
            time_vals,
            z,
            swd,
            levels = 50,
            colormap = :viridis,
            extendlow = :red,
        )
        graph_num += 1
    end
    Makie.Colorbar(
        figure[2, :],
        limits = (low_lim, up_lim),
        colormap = :viridis,
        label = "Shortwave Flux down [W/m²]",
        vertical = false,
    )
    Makie.display(figure)
    # Makie.save(figureName, figure)
end

# plots longwave flux contour plot over time (similar to last)
function plot_lw_flux_down_time(
    folder_paths,
    graph_titles,
    figureName,
    low_lim,
    up_lim,
)
    figure = Makie.Figure()
    Makie.empty!(figure)
    days_file_paths = create_day_file_paths(40, 200)
    time_vals = Float32[0.0, 40.0, 80.0, 120.0, 160.0, 200.0]
    _, diagnostics_example =
        read_hdf5_file(folder_paths[1] * days_file_paths[1])
    _, z... = vec(
        parent(
            Fields.coordinate_field(axes(diagnostics_example.lw_flux_down)).z,
        ),
    )
    plots_per_row = Int(length(folder_paths) / 2)
    graph_num = 1
    for folder_path in folder_paths
        time_count = 1
        lwd = Matrix{Float32}(undef, length(days_file_paths), length(z))
        for day in days_file_paths
            file_path = folder_path * day
            time, diagnostics = read_hdf5_file(file_path)
            _, lwd[time_count, :]... = vec(parent(diagnostics.lw_flux_down))
            time_count += 1
        end
        axis = Makie.Axis(
            figure[
                trunc(Int, (graph_num - 1) / plots_per_row),
                (graph_num - 1) % plots_per_row,
            ],
            xlabel = "Time [days]",
            ylabel = "Elevation [m]",
            title = graph_titles[graph_num],
            height = 160,
            width = 165 * 3 / plots_per_row,
            xlabelsize = 14.0,
            ylabelsize = 14.0,
            xminorticksize = 2.0,
            yminorticksize = 2.0,
            xticklabelsize = 10.0,
            yticklabelsize = 10.0,
        )
        Makie.contourf!(
            axis,
            time_vals,
            z,
            lwd,
            levels = 50,
            colormap = :viridis,
            extendlow = :red,
        )
        graph_num += 1
    end
    Makie.Colorbar(
        figure[2, :],
        limits = (low_lim, up_lim),
        colormap = :viridis,
        label = "Longwave Flux down [W/m²]",
        vertical = false,
    )
    Makie.display(figure)
    # Makie.save(figureName, figure)
end


# plots the surface temperature change over time of a single simulation
function plot_surface_temp_convergence(interval, max_day, folder_path)
    figure = Makie.Figure()
    days_file_paths = create_day_file_paths(interval, max_day)
    time_vals = create_time_vals(interval, max_day)
    time_count = 1
    sfc_temp_vec = zeros(length(time_vals))
    for day in days_file_paths
        file_path = folder_path * day
        time, diagnostics = read_hdf5_file(file_path)
        sfc_temp_vec[time_count] = vec(parent(diagnostics.sfc_temperature))[1]
        time_count += 1
    end
    axis = Makie.Axis(
        figure[1, 1],
        xlabel = "time [days]",
        ylabel = "surface temp [K]",
        title = "Surface Temp convergence with prognostic surface temperature",
    )
    Makie.lines!(axis, time_vals, sfc_temp_vec)
    Makie.display(figure)
    Makie.save("Surface_temp_convergence.png", figure)
end

# plots multiple temperature profiles at once (4)
function plot_four_temp_profiles(max_day, folder_path, subtitles, saveName, title)
    figure = Makie.Figure()
    position = 1
    axis = Makie.Axis(
        figure[1, 1],
        xlabel = "temperature [K]",
        ylabel = "height [m]",
        title = title,
    )
    for folder in folder_path
        file_path = folder * create_day_file_paths(max_day, max_day)[2]
        _, diagnostics = read_hdf5_file(file_path)
        sfc_t = vec(parent(diagnostics.sfc_temperature))[1]
        sfc_z = 0.0
        z = vec(
            parent(Fields.coordinate_field(axes(diagnostics.temperature)).z),
        )
        temp = vec(parent(diagnostics.temperature))
        if (position == 1)
            Makie.lines!(axis, temp, z, color = :green)
            Makie.scatter!(axis, sfc_t, sfc_z, color = :green)
        elseif (position == 2)
            Makie.lines!(axis, temp, z, color = :blue)
            Makie.scatter!(axis, sfc_t, sfc_z, color = :blue)
        elseif (position == 3)
            Makie.lines!(axis, temp, z, color = :orange)
            Makie.scatter!(axis, sfc_t, sfc_z, color = :orange)
        else
            Makie.lines!(axis, temp, z, color = :red)
            Makie.scatter!(axis, sfc_t, sfc_z, color = :red)
        end
        position += 1
    end
    Makie.Legend(
        figure[1, 2],
        [
            Makie.LineElement(color = :green, linestyle = nothing),
            Makie.LineElement(color = :blue, linestyle = nothing),
            Makie.LineElement(color = :orange, linestyle = nothing),
            Makie.LineElement(color = :red, linestyle = nothing),
            Makie.MarkerElement(
                color = :black,
                marker = :circle,
                markersize = 14,
            ),
        ],
        subtitles,
    )
    Makie.display(figure)
    Makie.save(saveName, figure)
end

# plots 6 temperature profiles at once
function plot_six_temp_profile(max_day, folder_path, subtitles, saveName, title)
    figure = Makie.Figure()
    position = 1
    axis = Makie.Axis(
        figure[1, 1],
        xlabel = "temperature [K]",
        ylabel = "height [m]",
        title = title,
    )
    for folder in folder_path
        file_path = folder * create_day_file_paths(max_day, max_day)[2]
        _, diagnostics = read_hdf5_file(file_path)
        sfc_t = vec(parent(diagnostics.sfc_temperature))[1]
        sfc_z = 0.0
        z = vec(
            parent(Fields.coordinate_field(axes(diagnostics.temperature)).z),
        )
        temp = vec(parent(diagnostics.temperature))
        if (position == 1)
            Makie.lines!(axis, temp, z, color = :green)
            Makie.scatter!(axis, sfc_t, sfc_z, color = :green)
        elseif (position == 2)
            Makie.lines!(axis, temp, z, color = :blue)
            Makie.scatter!(axis, sfc_t, sfc_z, color = :blue)
        elseif (position == 3)
            Makie.lines!(axis, temp, z, color = :black)
            Makie.scatter!(axis, sfc_t, sfc_z, color = :blue)
        elseif (position == 4)
            Makie.lines!(axis, temp, z, color = :yellow)
            Makie.scatter!(axis, sfc_t, sfc_z, color = :yellow)
        elseif (position == 5)
            Makie.lines!(axis, temp, z, color = :orange)
            Makie.scatter!(axis, sfc_t, sfc_z, color = :orange)
        else
            Makie.lines!(axis, temp, z, color = :red)
            Makie.scatter!(axis, sfc_t, sfc_z, color = :red)
        end
        position += 1
    end
    Makie.Legend(
        figure[1, 2],
        [
            Makie.LineElement(color = :green, linestyle = nothing),
            Makie.LineElement(color = :blue, linestyle = nothing),
            Makie.LineElement(color = :black, linestyle = nothing),
            Makie.LineElement(color = :yellow, linestyle = nothing),
            Makie.LineElement(color = :orange, linestyle = nothing),
            Makie.LineElement(color = :red, linestyle = nothing),
            Makie.MarkerElement(
                color = :black,
                marker = :circle,
                markersize = 14,
            ),
        ],
        subtitles,
    )
    Makie.display(figure)
    Makie.save(saveName, figure)
end

# plots how surface temperature converges for differing coefficients in the surface temperature equation
function plot_varied_convergence_surface_temp(folder_paths, interval, max_day)
    figure = Makie.Figure()
    days_file_paths = create_day_file_paths(interval, max_day)
    time_vals = create_time_vals(interval, max_day)
    axis = Makie.Axis(
        figure[1, 1],
        xlabel = "Day",
        ylabel = "Surface temperature [K]",
        title = "Convergence of Surface Temperature",
        subtitle = "At different mixed ocean layer depths",
    )
    legend_count = 1
    for folder in folder_paths
        time_count = 1
        sfc_temp_vec = zeros(length(time_vals))
        for day in days_file_paths
            file_path = folder * day
            _, diagnostics = read_hdf5_file(file_path)
            sfc_temp_vec[time_count] =
                vec(parent(diagnostics.sfc_temperature))[1]
            time_count += 1
        end
        Makie.lines!(axis, time_vals, sfc_temp_vec)
        legend_count += 1
    end
    # elem1 = [LineElement(color = :blue, linestyle = nothing), strokecolor = :black]
    Makie.Legend(
        figure[1, 2],
        [
            Makie.LineElement(color = :blue, linestyle = nothing),
            Makie.LineElement(color = :orange, linestyle = nothing),
            Makie.LineElement(color = :green, linestyle = nothing),
            Makie.LineElement(color = :red, linestyle = nothing),
        ],
        ["depth = 0.5 m", "depth = 2 m", "depth = 5 m", "depth = 20 m"],
    )
    Makie.display(figure)
    Makie.save("varied_surface_temp_convergence.png", figure)
end

# plots the temperature difference between two different simulations at a certain time
function plot_temp_profile_difference(path1, path2, subtitle, saveName)
    figure = Makie.Figure()
    Makie.empty!(figure)
    _, diagnostics1 = read_hdf5_file(path1)
    _, diagnostics2 = read_hdf5_file(path2)
    z = vec(parent(Fields.coordinate_field(axes(diagnostics1.temperature)).z))
    temp_diff =
        vec(parent(diagnostics1.temperature)) -
        vec(parent(diagnostics2.temperature))
    axis = Makie.Axis(
        figure[1, 1],
        xlabel = "Temperature [K]",
        ylabel = "Elevation [m]",
        title = "Difference Temperature between Prognostic and Prescribed Surface models",
        subtitle = subtitle,
    )
    Makie.lines!(axis, temp_diff, z)
    Makie.display(figure)
    Makie.save(saveName, figure)
end

# plots the temperature of different simulations, now with respect to pressure instead of height/altitude
function plot_temp_pressure_profile(file_path, saveName, title)
    figure = Makie.Figure()
    axis = Makie.Axis(
        figure[1, 1],
        xlabel = "temperature [K]",
        ylabel = "pressure [hPa]",
        title = title,
    )
    _, diagnostics = read_hdf5_file(file_path)
    sfc_t = vec(parent(diagnostics.sfc_temperature))[1]
    pressure = vec(parent(diagnostics.pressure))
    sfc_pressure = pressure[1]
    temp = vec(parent(diagnostics.temperature))
    Makie.lines!(axis, temp, pressure, color = :blue)
    Makie.scatter!(axis, sfc_t, sfc_pressure, color = :blue)
    Makie.display(figure)
    Makie.save(saveName, figure)
end

# plots the difference in temp between simulations (with respect to pressure)
function plot_temp_pressure_profile_difference(file_path1, file_path2, saveName, title)
    figure = Makie.Figure()
    axis = Makie.Axis(
        figure[1, 1],
        xlabel = "temperature difference [K]",
        ylabel = "height [m]", 
        title = title,
    )
    Makie.xlims!(-11.0, 13.0)
    _, diagnostics = read_hdf5_file(file_path1)
    _, diagnostics2 = read_hdf5_file(file_path2)
    sfc_t =
        vec(parent(diagnostics.sfc_temperature))[1] -
        vec(parent(diagnostics2.sfc_temperature))[1]
    z = vec(parent(Fields.coordinate_field(axes(diagnostics.temperature)).z))
    sfc_z = z[1]
    temp =
        vec(parent(diagnostics.temperature)) -
        vec(parent(diagnostics2.temperature))
    Makie.lines!(axis, temp, z, color = :red)
    Makie.scatter!(axis, sfc_t, sfc_z, color = :red)
    Makie.display(figure)
    Makie.save(saveName, figure)
end

# compares multiple temperature profiles from different simulations
function plot_temp_profile_differences(
    max_day,
    folder_path,
    subtitles,
    saveName,
    title,
    limit,
)
    figure = Makie.Figure()
    position = 1
    axis = Makie.Axis(
        figure[1, 1],
        xlabel = "temperature difference [K]",
        ylabel = "height [m]",
        title = title,
    )
    Makie.xlims!(-limit, limit)
    def_folder = folder_path[2]
    def_file_path = def_folder * create_day_file_paths(max_day, max_day)[2]
    _, def_diagnostics = read_hdf5_file(def_file_path)
    def_sfc_t = sfc_t = vec(parent(def_diagnostics.sfc_temperature))[1]
    def_temp = vec(parent(def_diagnostics.temperature))

    for folder in folder_path
        if (folder != def_folder)
            file_path = folder * create_day_file_paths(max_day, max_day)[2]
            _, diagnostics = read_hdf5_file(file_path)
            sfc_t = vec(parent(diagnostics.sfc_temperature))[1] - def_sfc_t
            sfc_z = 0.0
            z = vec(
                parent(
                    Fields.coordinate_field(axes(diagnostics.temperature)).z,
                ),
            )
            temp = vec(parent(diagnostics.temperature)) - def_temp
            if (position == 1)
                Makie.lines!(axis, temp, z, color = :green)
                Makie.scatter!(axis, sfc_t, sfc_z, color = :green)
            elseif (position == 2)
                Makie.lines!(axis, temp, z, color = :orange)
                Makie.scatter!(axis, sfc_t, sfc_z, color = :orange)
            else
                Makie.lines!(axis, temp, z, color = :red)
                Makie.scatter!(axis, sfc_t, sfc_z, color = :red)
            end
            position += 1
        end
    end
    Makie.Legend(
        figure[1, 2],
        [
            Makie.LineElement(color = :green, linestyle = nothing),
            Makie.LineElement(color = :orange, linestyle = nothing),
            Makie.LineElement(color = :red, linestyle = nothing),
            Makie.MarkerElement(
                color = :black,
                marker = :circle,
                markersize = 14,
            ),
        ],
        subtitles,
    )
    Makie.display(figure)
    Makie.save(saveName, figure)
end
