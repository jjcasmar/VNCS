#ifndef LOGGER_H
#define LOGGER_H

#include <spdlog/common.h>
#include <spdlog/spdlog.h>
#include <spdlog/logger.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/fmt/ostr.h>
#include <type_traits>

namespace VNCS
{
template <typename LoggerCategory>
struct StaticLogger {
    StaticLogger(const std::string &loggerName) {}

    static auto logger() { return spdlog::get(LoggerCategory::name); }

    template <typename FormatString, typename... Args>
    static void trace(const FormatString &fmt, Args &&...args)
    {
        const auto logger = spdlog::get(LoggerCategory::name);
        if (logger)
            logger->trace(fmt, std::forward<Args>(args)...);
    }

    template <typename FormatString, typename... Args>
    static void debug(const FormatString &fmt, Args &&...args)
    {
        const auto logger = spdlog::get(LoggerCategory::name);
        if (logger)
            logger->debug(fmt, std::forward<Args>(args)...);
    }

    template <typename FormatString, typename... Args>
    static void info(const FormatString &fmt, Args &&...args)
    {
        const auto logger = spdlog::get(LoggerCategory::name);
        if (logger)
            logger->info(fmt, std::forward<Args>(args)...);
    }

    template <typename FormatString, typename... Args>
    static void warn(const FormatString &fmt, Args &&...args)
    {
        const auto logger = spdlog::get(LoggerCategory::name);
        if (logger)
            logger->warn(fmt, std::forward<Args>(args)...);
    }

    template <typename FormatString, typename... Args>
    static void error(const FormatString &fmt, Args &&...args)
    {
        const auto logger = spdlog::get(LoggerCategory::name);
        if (logger)
            logger->error(fmt, std::forward<Args>(args)...);
    }

    template <typename FormatString, typename... Args>
    static void critical(const FormatString &fmt, Args &&...args)
    {
        const auto logger = spdlog::get(LoggerCategory::name);
        if (logger)
            logger->error(fmt, std::forward<Args>(args)...);
    }
};  // namespace Log

struct DynamicLogger {
    DynamicLogger(const std::string &loggerName, spdlog::level::level_enum level = spdlog::level::info)
        : m_loggerName(loggerName)
    {
        if (!spdlog::get(loggerName)) {
            auto logger = spdlog::stdout_color_mt(loggerName);
            logger->set_level(level);
        }
    }

    auto logger() { return spdlog::get(m_loggerName); }

    template <typename FormatString, typename... Args>
    void trace(const FormatString &fmt, Args &&...args)
    {
        const auto logger = spdlog::get(m_loggerName);
        if (logger)
            logger->trace(fmt, std::forward<Args>(args)...);
    }

    template <typename FormatString, typename... Args>
    void debug(const FormatString &fmt, Args &&...args)
    {
        const auto logger = spdlog::get(m_loggerName);
        if (logger)
            logger->debug(fmt, std::forward<Args>(args)...);
    }

    template <typename FormatString, typename... Args>
    void info(const FormatString &fmt, Args &&...args)
    {
        const auto logger = spdlog::get(m_loggerName);
        if (logger)
            logger->info(fmt, std::forward<Args>(args)...);
    }

    template <typename FormatString, typename... Args>
    void warn(const FormatString &fmt, Args &&...args)
    {
        const auto logger = spdlog::get(m_loggerName);
        if (logger)
            logger->warn(fmt, std::forward<Args>(args)...);
    }

    template <typename FormatString, typename... Args>
    void error(const FormatString &fmt, Args &&...args)
    {
        const auto logger = spdlog::get(m_loggerName);
        if (logger)
            logger->error(fmt, std::forward<Args>(args)...);
    }

    template <typename FormatString, typename... Args>
    void critical(const FormatString &fmt, Args &&...args)
    {
        const auto logger = spdlog::get(m_loggerName);
        if (logger)
            logger->error(fmt, std::forward<Args>(args)...);
    }

private:
    std::string m_loggerName;
};  // namespace Log

template <typename T>
struct LoggerImplementation {
    using type = StaticLogger<T>;
};

template <>
struct LoggerImplementation<DynamicLogger> {
    using type = DynamicLogger;
};

template <typename LoggerCategory = DynamicLogger>
class Logger : public LoggerImplementation<LoggerCategory>::type
{
public:
    template <typename... Args>
    Logger(Args &&...args)
        : LoggerImplementation<LoggerCategory>::type(std::forward<Args>(args)...)
    {
    }
};
}  // namespace VNCS

#endif  // LOGGER_H
