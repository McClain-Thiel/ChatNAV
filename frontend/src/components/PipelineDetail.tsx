import type { ModuleState } from '../api/types';
import { formatDuration } from '../api/helpers';
import FileViewer from './FileViewer';

interface Props {
  module: ModuleState;
  totalModules: number;
  jobId: string;
}

export default function PipelineDetail({ module, totalModules, jobId }: Props) {
  const isRunning = module.status === 'running';
  const isComplete = module.status === 'complete';

  return (
    <div className="animate-reveal" style={{ maxWidth: 640, margin: '0 auto', padding: '44px 52px' }}>
      {/* Running status line */}
      {isRunning && module.progress_pct != null && (
        <div style={{
          display: 'flex',
          alignItems: 'center',
          gap: 8,
          marginBottom: 20,
        }}>
          <span className="animate-blink" style={{
            fontFamily: "'Inconsolata', monospace",
            fontSize: 9,
            letterSpacing: 2,
            textTransform: 'uppercase',
            color: 'var(--status-running)',
          }}>
            ● Running — {module.progress_pct}% complete
          </span>
          <div style={{ flex: 1, height: 1, background: 'var(--border-light)' }} />
        </div>
      )}

      {/* Step kicker */}
      <div style={{
        fontFamily: "'Inconsolata', monospace",
        fontSize: 8.5,
        color: 'var(--text-muted)',
        letterSpacing: 3,
        textTransform: 'uppercase',
        marginBottom: 8,
      }}>
        Step {String(module.number).padStart(2, '0')} of {totalModules}
      </div>

      {/* Headline */}
      <h1 style={{
        fontFamily: "'EB Garamond', serif",
        fontSize: 34,
        fontWeight: 500,
        color: 'var(--text-primary)',
        marginBottom: 6,
        lineHeight: 1.2,
      }}>
        {module.name}
      </h1>

      {/* Deck */}
      <p style={{
        fontFamily: "'EB Garamond', serif",
        fontStyle: 'italic',
        fontSize: 17,
        color: 'var(--text-meta)',
        marginBottom: 20,
        lineHeight: 1.4,
      }}>
        {module.subtitle}
      </p>

      {/* Byline rule */}
      <div style={{ marginBottom: 28 }}>
        <div style={{ borderTop: '2px solid var(--border-heavy)' }} />
        <div style={{
          display: 'flex',
          gap: 24,
          padding: '6px 0',
          fontFamily: "'Inconsolata', monospace",
          fontSize: 9,
          color: 'var(--text-muted)',
        }}>
          <span>Tool: {module.tool}</span>
          {isComplete && module.duration_seconds != null && (
            <span style={{ color: 'var(--status-complete)' }}>
              Completed in {formatDuration(module.duration_seconds)}
            </span>
          )}
        </div>
        <div style={{ borderTop: '1px solid var(--border-mid)' }} />
      </div>

      {/* Progress bar (running only) */}
      {isRunning && module.progress_pct != null && (
        <div style={{ marginBottom: 28 }}>
          <div style={{
            height: 3,
            background: 'var(--border-light)',
            borderRadius: 1,
          }}>
            <div style={{
              height: '100%',
              width: `${module.progress_pct}%`,
              background: 'var(--status-running)',
              borderRadius: 1,
              transition: 'width 0.5s ease',
            }} />
          </div>
          {module.progress_detail && (
            <div style={{
              fontFamily: "'Inconsolata', monospace",
              fontSize: 9,
              color: 'var(--text-faint)',
              marginTop: 6,
            }}>
              {module.progress_detail}
            </div>
          )}
        </div>
      )}

      {/* Body copy */}
      {module.description && (
        <div style={{
          fontFamily: "'Spectral', serif",
          fontSize: 15.5,
          lineHeight: 1.9,
          color: 'var(--text-secondary)',
          textAlign: 'justify',
          marginBottom: 28,
        }}>
          {module.description}
        </div>
      )}

      {/* Scientific significance */}
      {module.significance && (
        <div style={{
          borderLeft: '2px solid var(--border-heavy)',
          paddingLeft: 18,
          marginBottom: 32,
        }}>
          <div style={{
            fontFamily: "'Inconsolata', monospace",
            fontSize: 8.5,
            letterSpacing: 2.5,
            textTransform: 'uppercase',
            color: 'var(--text-meta)',
            marginBottom: 8,
          }}>
            Scientific Significance
          </div>
          <div style={{
            fontFamily: "'Spectral', serif",
            fontStyle: 'italic',
            fontSize: 15,
            lineHeight: 1.85,
            color: 'var(--text-secondary)',
          }}>
            {module.significance}
          </div>
        </div>
      )}

      {/* Output files */}
      {module.outputs.length > 0 && (
        <div style={{
          borderTop: '1px solid var(--border-mid)',
          paddingTop: 22,
          marginBottom: 24,
        }}>
          <div style={{
            fontFamily: "'Inconsolata', monospace",
            fontSize: 8.5,
            letterSpacing: 2.5,
            textTransform: 'uppercase',
            color: 'var(--text-muted)',
            marginBottom: 12,
          }}>
            Output Files
          </div>
          {isComplete && module.output_dir ? (
            module.outputs.map(file => (
              <FileViewer
                key={file}
                jobId={jobId}
                filePath={file}
                outputDir={module.output_dir!}
              />
            ))
          ) : (
            <div style={{
              display: 'grid',
              gridTemplateColumns: '1fr 1fr',
              gap: '6px 16px',
            }}>
              {module.outputs.map(file => (
                <div key={file} style={{
                  fontFamily: "'Inconsolata', monospace",
                  fontSize: 10,
                  color: 'var(--text-secondary)',
                  display: 'flex',
                  alignItems: 'center',
                  gap: 6,
                }}>
                  <span style={{ color: isComplete ? 'var(--status-complete)' : 'var(--text-dim)' }}>
                    {isComplete ? '✓' : '○'}
                  </span>
                  {file}
                </div>
              ))}
            </div>
          )}
        </div>
      )}

      {/* References */}
      {module.references.length > 0 && (
        <div style={{
          borderTop: '1px solid var(--border-mid)',
          paddingTop: 20,
        }}>
          <div style={{
            fontFamily: "'Inconsolata', monospace",
            fontSize: 8.5,
            letterSpacing: 2.5,
            textTransform: 'uppercase',
            color: 'var(--text-muted)',
            marginBottom: 10,
          }}>
            References
          </div>
          {module.references.map((ref, i) => (
            <div key={i} style={{
              fontFamily: "'Spectral', serif",
              fontStyle: 'italic',
              fontSize: 11,
              color: 'var(--text-meta)',
              marginBottom: 4,
            }}>
              <span style={{ fontStyle: 'normal', color: 'var(--text-muted)', marginRight: 4 }}>
                [{i + 1}]
              </span>
              {ref}
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
