# Slack 플러그인

이 repository에는 Slack을 Cursor IDE 및 Claude Code와 통합하는 데 필요한 configuration이 포함되어 있습니다. 이 plugin을 사용하면 agent가 Slack workspace와 직접 상호작용해 message search, communication send, canvas management 등을 natural language로 수행할 수 있습니다.

## 기능

Slack MCP server는 다음 기능을 제공합니다:

- **Search**: Public/private message, file, user, channel을 찾습니다
- **Messaging**: Message를 보내고 channel history와 thread conversation을 가져옵니다
- **Canvas**: Formatted document를 만들고 공유하며 markdown으로 export합니다
- **User Management**: Custom field와 status information을 포함한 user profile을 가져옵니다

## 사전 요구 사항

Slack MCP server를 설정하기 전에 다음을 확인하세요.

- Cursor IDE 또는 Claude Code CLI가 설치되어 있어야 합니다
- Workspace admin이 MCP integration을 승인한 Slack workspace에 접근할 수 있어야 합니다

## 설치

IDE에 맞는 설치 방법을 선택하세요.

### Claude Code

Claude Code CLI를 사용한다면 local clone으로 plugin을 설치할 수 있습니다.

```bash
git clone https://github.com/slackapi/slack-mcp-plugin.git
cd slack-mcp-plugin
claude --plugin-dir ./
```

Plugin이 load되면 Slack MCP server가 자동으로 설정됩니다. OAuth로 Slack workspace에 authenticate하라는 prompt가 표시됩니다.

Claude plugin은 다음 MCP configuration(`.mcp.json`)을 사용합니다.

```json
{
  "mcpServers": {
    "slack": {
      "type": "http",
      "url": "https://mcp.slack.com/mcp",
      "oauth": {
        "clientId": "1601185624273.8899143856786",
        "callbackPort": 3118
      }
    }
  }
}
```

### Cursor

아래 Add to Cursor button을 사용하거나, 다음 단계에 따라 Cursor에서 Slack MCP server를 수동으로 설정할 수 있습니다.

[![Install MCP Server](https://cursor.com/deeplink/mcp-install-dark.svg)](https://cursor.com/en-US/install-mcp?name=slack&config=eyJ1cmwiOiJodHRwczovL21jcC5zbGFjay5jb20vbWNwIiwiYXV0aCI6eyJDTElFTlRfSUQiOiIzNjYwNzUzMTkyNjI2Ljg5MDM0NjkyMjg5ODIifX0%3D)

#### Step 1: Cursor settings 열기

**Cursor → Settings → Cursor Settings**로 이동합니다. macOS에서는 `Cmd+,`, Windows/Linux에서는 `Ctrl+,` shortcut을 사용할 수 있습니다.

#### Step 2: MCP tab으로 이동

Settings interface에서 **MCP** tab을 클릭해 MCP server configuration에 접근합니다.

#### Step 3: Slack MCP configuration 추가

Remote Slack MCP server에 연결하려면 다음 configuration을 추가합니다.

```json
{
  "mcpServers": {
    "slack": {
      "url": "https://mcp.slack.com/mcp",
      "auth": {
        "CLIENT_ID": "3660753192626.8903469228982"
      }
    }
  }
}
```

Configuration을 저장합니다. 추가 후 connect button이 보이면 클릭해 Slack Workspace에 authenticate합니다.

## 사용 예시

설정이 끝나면 자연어로 AI assistant를 통해 Slack과 상호작용할 수 있습니다.

- **Search messages**: "지난주 product launch 관련 message를 찾아줘"
- **Send messages**: "#general channel에 deployment가 완료됐다고 보내줘"
- **Find users**: "john@example.com email을 가진 user가 누구야?"
- **Access threads**: "그 message의 conversation thread를 보여줘"
- **Create canvases**: "Meeting note로 canvas document를 만들어줘"

## 문서 및 resource

- [Official Slack MCP Server Documentation](https://docs.slack.dev/ai/mcp-server/)

## 참고 및 제한 사항

- **Remote server only**: 이 configuration은 Slack hosted MCP server에 연결합니다. Local installation은 필요하지 않으며 지원되지 않습니다.
- **Admin approval required**: 이 기능을 사용하려면 Slack workspace administrator가 MCP integration을 승인해야 합니다.

## 질문 또는 issue

Slack MCP server 질문이나 integration issue는 [official Slack documentation](https://docs.slack.dev/ai/mcp-server/)을 참고하거나 workspace administrator에게 문의하세요.
