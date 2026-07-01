//
// Created by Niharika Gauraha on 2026-07-01.
//

#ifndef MERGELLL_THREADPOOL_H
#define MERGELLL_THREADPOOL_H
#include <future>
#include <type_traits>
#include <functional>

class ThreadPool {
public:
    ThreadPool(int threads);
    ~ThreadPool();

    template<class F>
    auto enqueue(F&& f) -> std::future<std::invoke_result_t<F>>;

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex m;
    std::condition_variable cv;
    bool stop;
};
template <class F>
auto ThreadPool::enqueue(F&& f)
    -> std::future<std::invoke_result_t<F>>
{
    using ReturnType = std::invoke_result_t<F>;

    auto task = std::make_shared<std::packaged_task<ReturnType()>>(
        std::forward<F>(f));

    std::future<ReturnType> future = task->get_future();

    {
        std::lock_guard<std::mutex> lock(m);
        tasks.emplace([task]() {
            (*task)();
        });
    }

    cv.notify_one();
    return future;
}

ThreadPool::~ThreadPool() {
    stop = true;
    cv.notify_all();

    for (auto &t : workers)
        if (t.joinable())
            t.join();
}

ThreadPool::ThreadPool(int threads) : stop(false) {
    for (int i = 0; i < threads; ++i) {
        workers.emplace_back([this] {
            while (true) {
                std::function<void()> task;

                {
                    std::unique_lock<std::mutex> lock(m);

                    cv.wait(lock, [this] {
                        return stop || !tasks.empty();
                    });

                    if (stop && tasks.empty())
                        return;

                    task = std::move(tasks.front());
                    tasks.pop();
                }

                task();
            }
        });
    }
}

#endif //MERGELLL_THREADPOOL_H
