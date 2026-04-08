import { useState, useEffect, useCallback, useRef } from 'react';
import type { JobDetail } from '../api/types';
import { getJob } from '../api/client';

export function useJob(jobId: string | undefined, pollInterval = 15_000) {
  const [job, setJob] = useState<JobDetail | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const jobRef = useRef<JobDetail | null>(null);

  const refresh = useCallback(async () => {
    if (!jobId) return;
    try {
      const data = await getJob(jobId);
      setJob(data);
      jobRef.current = data;
      setError(null);
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to load job');
    } finally {
      setLoading(false);
    }
  }, [jobId]);

  useEffect(() => {
    setLoading(true);
    setJob(null);
    jobRef.current = null;
    refresh();

    const id = setInterval(() => {
      const status = jobRef.current?.status;
      if (status === 'completed' || status === 'failed') return;
      refresh();
    }, pollInterval);

    return () => clearInterval(id);
  }, [jobId, refresh, pollInterval]);

  return { job, loading, error, refresh };
}
