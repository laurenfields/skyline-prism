using System;
using System.Collections.Generic;
using System.Net.Http;
using System.Text;
using System.Text.Json;

namespace SkylinePrism.Core.DifferentialAnalysis.Enrichment;

/// <summary>
/// Default <see cref="IJsonPoster"/>: POSTs the payload as JSON over HTTP and returns the parsed
/// response. Used by the enrichment clients for their live g:Profiler / Open Targets calls; tests
/// inject a fake instead. Reuses a single <see cref="HttpClient"/> (safe to share).
/// </summary>
public sealed class HttpJsonPoster : IJsonPoster, IDisposable
{
    private readonly HttpClient _client;
    private readonly bool _ownsClient;

    /// <summary>
    /// The default request timeout. Two minutes rather than the 30 seconds this started at, because
    /// the g:Profiler query now asks for per-gene evidence (<c>no_evidences = false</c>, which is what
    /// lets a term name its member proteins) and the response carries one evidence list per query gene
    /// per term. A proteome-scale run - a thousand significant genes against every quantified gene,
    /// over five sources - returns enough that the original timeout was sized for a response the code
    /// no longer asks for. This is a one-off interactive query, never a hot path, so waiting is
    /// cheaper than failing.
    /// </summary>
    public static readonly TimeSpan DefaultTimeout = TimeSpan.FromMinutes(2);

    /// <summary>Create a poster with its own <see cref="HttpClient"/> and the given request timeout.</summary>
    public HttpJsonPoster(TimeSpan? timeout = null)
    {
        _client = new HttpClient { Timeout = timeout ?? DefaultTimeout };
        _ownsClient = true;
    }

    /// <summary>Create a poster over an existing <see cref="HttpClient"/> (not disposed by this poster).</summary>
    public HttpJsonPoster(HttpClient client)
    {
        _client = client;
        _ownsClient = false;
    }

    /// <inheritdoc/>
    public JsonElement Post(string url, IReadOnlyDictionary<string, object?> payload)
    {
        var json = JsonSerializer.Serialize(payload);
        using var content = new StringContent(json, Encoding.UTF8, "application/json");
        using var response = _client.PostAsync(url, content).GetAwaiter().GetResult();
        response.EnsureSuccessStatusCode();
        var body = response.Content.ReadAsStringAsync().GetAwaiter().GetResult();
        using var doc = JsonDocument.Parse(body);
        return doc.RootElement.Clone();
    }

    /// <inheritdoc/>
    public void Dispose()
    {
        if (_ownsClient)
            _client.Dispose();
    }
}
